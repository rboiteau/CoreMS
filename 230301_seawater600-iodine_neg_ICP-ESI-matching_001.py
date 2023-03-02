import re
import os 
from tempfile import tempdir
import time
from turtle import color
import numpy as np
import warnings
from datetime import date, datetime
from pathlib import Path
import sys
import tqdm
import pandas as pd
import matplotlib.pyplot as plt
import pickle
sys.path.append("./")
warnings.filterwarnings("ignore")

from corems.mass_spectra.input import rawFileReader
from corems.molecular_id.factory.classification import HeteroatomsClassification, Labels
from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.encapsulation.factory.parameters import MSParameters
from corems.encapsulation.constant import Atoms
from corems.mass_spectrum.calc.Calibration import MzDomainCalibration
import corems.lc_icpms_ftms.calc.lc_icrms_helpers as lcmsfns

"""
Author:     Christian Dewey
Run on:     3 Mar 23

CoreMS run script for correlating ICPMS, ESIMS data

    Dataset:        wastewater_SE_RC_ENV

    Instruments:    Keck iCAP RQ, Keck IQ-X

    Analysis date:  28 Feb 23 (IQ-X); 23 Feb 23 (iCAP) 

    Analysis notes: attempting to identify I features 
"""

def plotICPMS(icp_data_file = None, elements = ['115In'], offset = 0, ax = None):
    # offset in seconds, required shift of ICPMS data to align with ESIMS data, based on b12 peak 
    if ax == None:
        _, ax = plt.subplots()

    
    icp_data = pd.read_csv(icp_data_file)
    try:
        etime = 'Time ' + elements[0]
        icp_subset = icp_data[[elements[0],etime]]
    except:
        icp_data = pd.read_csv(icp_data_file, sep=';', skiprows=1)
        
    icp_data.dropna(inplace=True)

    for element in elements:
        etime = 'Time ' + element
        icp_subset = icp_data[[element,etime]]
        icp_subset[etime] = (icp_subset[etime]+ offset) / 60
        ax.plot(icp_subset[etime], icp_subset[element]/1e4, label=element)
    ax.set_ylabel('ICP-MS Intensity (10$^4$ cps)')
    ax.set_xlabel('Time (min)')
    ax.set_title(icp_data_file.split('/')[-1])
    ax.set_ylim(bottom = 0)
    ax.legend(frameon = False)
    return ax

def mouse_event(event):
    print('x: {} and y: {}'.format(event.xdata, event.ydata))
    ix, iy = event.xdata, event.ydata
    coords.append((ix,iy))
    if len(coords) == 2:
        fig.canvas.mpl_disconnect(cid)
        plt.close('all')


def subset_icpdata(icp_data_file = None, heteroAtom = '127I', timerange = [0,1], offset = 0):

    icp_data = pd.read_csv(icp_data_file)
    try:
        etime = 'Time ' + heteroAtom
        icp_subset = icp_data[[heteroAtom,etime]]
    except:
        icp_data = pd.read_csv(icp_data_file, sep=';', skiprows=1)
    icp_data.dropna(inplace=True)

    timestart = timerange[0]
    timestop = timerange[1]
    element = heteroAtom
    etime = 'Time ' + heteroAtom
    icp_subset = icp_data[[element,etime]]
    icp_subset[etime] = (icp_subset[etime] + offset) / 60 
    icp_subset = icp_subset[icp_subset[etime].between(timestart,timestop)]

    return icp_subset


def get_eics(esi_parser=None,assignments = None, timerange = None):

    timestart = timerange[0]
    timestop = timerange[1]
    tic = esi_parser.get_tic(ms_type = 'MS')[0]
    tic_df=pd.DataFrame({'time': tic.time,'scan': tic.scans})
    scans=tic_df[tic_df.time.between(timestart,timestop)].scan.tolist()
    AverageMS = esi_parser.get_average_mass_spectrum_by_scanlist(scans)

    EICdic = {}
    pbar = tqdm.tqdm(assignments['m/z'], desc="Getting EICs")
    
    for mz in pbar:   
        #   print('AverageMS mz:' + str(mz))
        EIC=esi_parser.get_eics(target_mzs=[mz],tic_data={},peak_detection=False,smooth=False)
        EICdic[mz]=EIC[0][mz]

    return EICdic, AverageMS


def interpolate(esi_parser = None, icpsub=None, heteroAtom = '127I', timerange = [0,1]):
    ###Interpolate LC-ICPMS data to obtain times matching ESI data

    timestart = timerange[0]
    timestop = timerange[1]
    etime = 'Time ' + heteroAtom
    tic = esi_parser.get_tic(ms_type = 'MS')[0]
    tic_df=pd.DataFrame({'time': tic.time,'scan': tic.scans})
    
    times=tic_df[tic_df.time.between(timestart,timestop)].time.tolist()

    icpsubset2 = icpsub

    pbar = tqdm.tqdm(times, desc="Subsetting ICPMS data" )
    
    for i in pbar:
        icpsubset2.loc[-1]=['NaN',i]
        icpsubset2 = icpsubset2.sort_index().reset_index(drop=True)

    
    icpsubset2=icpsubset2.sort_values(by=etime)
    icpsubset2=icpsubset2.astype(float)
    icpsubset3=icpsubset2.interpolate()
    
    icp_interp=pd.DataFrame()
    pbar = tqdm.tqdm(times,desc="Interpolating ICPMS data")
    
    for i in pbar:
        icp_interp=icp_interp.append(icpsubset3[icpsubset3[etime]==i])
        
    return icp_interp


def correlate(icp_interp=None,EICdic=None,heteroAtom='127I',assignments = None, timerange=[0,1],threshold = 0.5):
    element = heteroAtom
    etime = 'Time ' + heteroAtom
    timestart = timerange[0]
    timestop = timerange[1]

    mscorr={}

    EICcurr=pd.DataFrame(index=icp_interp[etime],columns=['EIC',element])
    EICcurr[element]=icp_interp[element].array

    pbar = tqdm.tqdm(EICdic.keys(),desc="Running correlation")
    
    for mz in pbar:
        EIC=pd.DataFrame({'EIC':EICdic[mz].eic,'Time':EICdic[mz].time})
        EIC_sub=EIC[EIC['Time'].between(timestart,timestop)]

        EICcurr['EIC']=EIC_sub['EIC'].values

        corvalue=EICcurr.corr(method='pearson')
        mscorr[mz]=corvalue.EIC[element]**2
        
    mzs_corr = pd.DataFrame.from_dict(mscorr,orient='index',columns=['corr'])

    assignments=assignments.sort_values(by=['m/z'])

    holder = pd.DataFrame( index = range(len(assignments['m/z'])), columns = ['m/z', 'corr'])  #index = assignments['Index'],

    holder = np.zeros((len(assignments['m/z']),2))

    for mz,row in zip(assignments['m/z'], range(len(assignments['m/z']))):

        holder[row,1] = mzs_corr[mzs_corr.index == mz].iloc[0]
        holder[row,0] = mz

    pdholder = pd.DataFrame(holder, columns = ['m/z', 'corr'])

    assignments.insert(4,'corr',pdholder['corr'].values)

    assignments.insert(5,'m/z-2',pdholder['m/z'].values)

    match = re.findall(r'[A-Za-z]+|\d+', heteroAtom)
    heteroAtom = match[1]
    
    results=assignments[assignments['corr']>threshold].filter(['m/z','Calibrated m/z','corr','Peak Height','Confidence Score','m/z Error (ppm)','Molecular Formula',heteroAtom])

    bestresults=results[(results[heteroAtom]>=1)]

    return results, bestresults


def plot_EIC_ICPMS(eic,mz,icp_data,element, trange, ax = None):

    if ax == None:
        _, ax = plt.subplots()

    icp_data.dropna(inplace=True)

    etime = 'Time ' + element
    icp_subset = icp_data[[element,etime]]
    icp_subset[etime] = (icp_subset[etime])
    ax.plot(icp_subset[etime], icp_subset[element]/max(icp_subset[element]), label=element)

    ax.plot(eic['Time'], eic['EIC'] / max(eic['EIC']), label='%.4f'%mz)
      
    ax.set_ylabel('Normalized Intensity')
    ax.set_xlabel('Time (min)')
    ax.set_ylim(bottom = 0)
    ax.legend(frameon = False)
    ax.set_xlim(trange[0],trange[1])
    return ax
    

def assign_formula(esiparser, trange, refmasslist=None,cal_ppm_threshold=(-1,1),charge=1):
        
    corder=2
    timestart=trange[0]
    timestop=trange[1]

    parser = esiparser

    tic=parser.get_tic(ms_type='MS')[0]
    tic_df=pd.DataFrame({'time': tic.time,'scan': tic.scans})
    
    scans=tic_df[tic_df.time.between(timestart,timestop)].scan.tolist()

    mass_spectrum = parser.get_average_mass_spectrum_by_scanlist(scans)    
    mass_spectrum.molecular_search_settings.ion_charge = charge

    if refmasslist:
        mass_spectrum.settings.min_calib_ppm_error = 10
        mass_spectrum.settings.max_calib_ppm_error = -10
        calfn = MzDomainCalibration(mass_spectrum, refmasslist)
        ref_mass_list_fmt = calfn.load_ref_mass_list(refmasslist)

        imzmeas, mzrefs = calfn.find_calibration_points(mass_spectrum, ref_mass_list_fmt,
                                                    calib_ppm_error_threshold=cal_ppm_threshold,
                                                    calib_snr_threshold=3)

        calfn.recalibrate_mass_spectrum(mass_spectrum, imzmeas, mzrefs, order=corder)


    SearchMolecularFormulas(mass_spectrum, first_hit=False).run_worker_mass_spectrum()

    mass_spectrum.percentile_assigned(report_error=True)

    assignments=mass_spectrum.to_dataframe()

    assignments['Time']=timestart

    return(assignments)   


def setAssingmentParams():
    # set assignment parameters
    MSParameters.mass_spectrum.threshold_method = 'signal_noise'
    MSParameters.mass_spectrum.s2n_threshold = 3
    MSParameters.ms_peak.peak_min_prominence_percent = 0.001

    MSParameters.molecular_search.error_method = 'None'
    MSParameters.molecular_search.min_ppm_error = -5.0
    MSParameters.molecular_search.max_ppm_error = 5.0

    MSParameters.molecular_search.isProtonated = True
    MSParameters.molecular_search.isRadical = False
    MSParameters.molecular_search.isAdduct = False

    MSParameters.molecular_search.score_method = "prob_score"
    MSParameters.molecular_search.output_score_method = "prob_score"

    MSParameters.molecular_search.url_database = 'postgresql+psycopg2://coremsappdb:coremsapppnnl@localhost:5432/coremsapp'
    MSParameters.molecular_search.min_dbe = -1
    MSParameters.molecular_search.max_dbe = 20

    MSParameters.molecular_search.usedAtoms['C'] = (1,50)
    MSParameters.molecular_search.usedAtoms['H'] = (4,100)
    MSParameters.molecular_search.usedAtoms['O'] = (1,20)
    MSParameters.molecular_search.usedAtoms['N'] = (0,4)
    MSParameters.molecular_search.usedAtoms['S'] = (0,4)
    MSParameters.molecular_search.usedAtoms['P'] = (0,4)
    MSParameters.molecular_search.usedAtoms['Cl'] = (0,2)
    MSParameters.molecular_search.usedAtoms['I'] = (0,1)
    MSParameters.molecular_search.usedAtoms['Br'] = (0,1)

    MSParameters.molecular_search.used_atom_valences={'C': (4),
                     '13C': (4),
                     'N': (3),
                     'O': (2),
                     'H': (1),
                     'Cl': (1, 0),
                     'I':(5),
                     'S': (2),
                     'P': (3, 5, 4, 2, 1),
                     'Br': (1, 0)}

if __name__ == '__main__':
    startdt = datetime.now()
    start = time.time()  #for duration
   
    drive_dir = '/Volumes/Samsung_T5/ESI-ICP-MS/seawater-iodine/'
    svdir=drive_dir
    
    mzref = "/Users/christiandewey/CoreMS/db/Hawkes_neg.ref"
    
    heteroAtom = '127I'

    offset = -38 # ICPMS time (s) + offset (s) = ESI time (s)




    esifile_name = '220822_CTD27_600m2.raw'
    data_dir = drive_dir
    esifile = data_dir+esifile_name
    esiparser = rawFileReader.ImportMassSpectraThermoMSFileReader(esifile)
    icpmsfile = data_dir + 'CTD27_600m.csv'
    results_fname='230301_CWD_220822_CTD27_600m2_001.csv'

      
    '''### plot offset ICPMS data and select peak for matching 
                fig, ax = plt.subplots()
                global coords
                coords = []
                cid = fig.canvas.mpl_connect('button_press_event', mouse_event)
                ax = plotICPMS(icpmsfile,['127I', '59Co'], offset, ax)
                ax.set_xlim(0,10)
                plt.show()
            
                
                ### subset ICP data 
                trange = [coords[0][0], coords[1][0]]'''

    trange = [44,47]  



    icpsub = subset_icpdata(icp_data_file=icpmsfile, heteroAtom='127I', timerange=trange, offset = offset)
  
    
    ### interpolate ICP data to match ESI data time points
    interpolated_ICP = interpolate(esi_parser=esiparser,icpsub=icpsub,heteroAtom='127I',timerange=trange)

    
    ### run formula assignment 
    setAssingmentParams()
    print(os.getcwd())
    assignments = assign_formula(esiparser,trange,mzref,cal_ppm_threshold=(-10,10),charge=1)
    assignments.to_csv(svdir+results_fname)
    all_results = pd.read_csv(svdir+results_fname)


    ### get EICS
    EICs, avMS = get_eics(esi_parser=esiparser,assignments = assignments, timerange=trange)



    with open(svdir + 'eics_'+results_fname + '.p', 'wb') as fp:
                    pickle.dump(EICs, fp, protocol=pickle.HIGHEST_PROTOCOL)
            
    with open(svdir + 'eics_'+results_fname + '.p', 'rb') as fp:
        EICs = pickle.load(fp)
    ### run correlation and filter assignment results
    corr_assignments, corr_assignments_with_hetero = correlate(icp_interp=interpolated_ICP,EICdic=EICs,heteroAtom='127I',assignments=assignments,timerange=trange,threshold = 0.2) 

    corr_assignments.to_csv(svdir+'corr_'+results_fname)
    
    for mz in corr_assignments['m/z']:
        fig, ax = plt.subplots()
        eic_df=pd.DataFrame({'EIC':EICs[mz].eic,'Time':EICs[mz].time})
        eic_sub=eic_df[eic_df['Time'].between(trange[0],trange[1])]
        #def plot_EIC_ICPMS(eic,mz,icp_data,element, trange, ax = None):
        plot_EIC_ICPMS(eic = eic_sub,mz = mz,icp_data = interpolated_ICP, element='127I', trange=trange, ax=ax)
        mf = corr_assignments[corr_assignments['m/z'] == mz]['Molecular Formula'].iloc[0]

        if pd.isna(mf):
            mf = 'Unassigned'
            mz_err = ''
        else:
            mz_err = 'm/z Error (ppm): %.3f' %(corr_assignments[corr_assignments['m/z'] == mz]['m/z Error (ppm)'].iloc[0])
        ax.text(0.05,0.95, mf, transform=ax.transAxes)
        ax.text(0.05,0.90, mz_err, transform=ax.transAxes)
        plt.savefig(svdir + results_fname.split('.')[0] + '_' + mf + '_' +'%.4f.pdf' %(mz))

    elapsed_time_sec = (time.time() - start) 

    if elapsed_time_sec > 3600:
        elapsed_time = elapsed_time_sec/3600
        unit = 'hr'
    else:
        elapsed_time = elapsed_time_sec / 60
        unit = 'min'

    enddt = datetime.now()
    startdt_str = startdt.strftime("%d-%b-%Y %H:%M:%S")
    enddt_str = enddt.strftime("%d-%b-%Y %H:%M:%S")

    print('\nAssignment took %.2f %s to complete' %(elapsed_time, unit))
    print('Started ' + startdt_str )
    print('Finished ' + enddt_str )














