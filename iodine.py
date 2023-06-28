import re
import os 
import time
import numpy as np
import warnings
from datetime import datetime
import sys
import tqdm
import pandas as pd
import matplotlib.pyplot as plt
import pickle

sys.path.append("./")
warnings.filterwarnings("ignore")

from corems.mass_spectra.input import rawFileReader
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.encapsulation.factory.parameters import MSParameters
from corems.mass_spectrum.calc.Calibration import MzDomainCalibration

"""
Author:     Christian Dewey
Run on:     26 June 23

CoreMS run script for correlating ICPMS, ESIMS data
"""



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


def get_eics(esi_file=None, assignments = None, timerange = None):
    esi_parser = rawFileReader.ImportMassSpectraThermoMSFileReader(esi_file)
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


def interpolate(esi_file = None, icpsub=None, heteroAtom = '127I', timerange = [0,1]):
    ###Interpolate LC-ICPMS data to obtain times matching ESI data
    esi_parser = rawFileReader.ImportMassSpectraThermoMSFileReader(esi_file)
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
    

def assign_formula(file, times, interval, refmasslist): 
    MSParameters.molecular_search.error_method = 'None'
    MSParameters.molecular_search.min_ppm_error = -5
    MSParameters.molecular_search.max_ppm_error = 5

    MSParameters.molecular_search.min_ion_charge = 1
    MSParameters.molecular_search.max_ion_charge = 2
    MSParameters.molecular_search.default_ion_charge = 1

    MSParameters.molecular_search.db_chunk_size = 500

    MSParameters.mass_spectrum.min_calib_ppm_error = -8
    MSParameters.mass_spectrum.max_calib_ppm_error = 0
    MSParameters.mass_spectrum.calib_pol_order = 2
    
    MSParameters.mass_spectrum.min_picking_mz=100
    MSParameters.mass_spectrum.max_picking_mz=900
    MSParameters.mass_spectrum.threshold_method = 'log'
    MSParameters.mass_spectrum.log_nsigma=400

    MSParameters.ms_peak.peak_min_prominence_percent = 0.001

    MSParameters.molecular_search.url_database = 'postgresql+psycopg2://coremsappdb:coremsapppnnl@localhost:5432/coremsapp' #'postgresql+psycopg2://coremsappdb:coremsapppnnl@molformdb-1:5432/coremsapp'
    MSParameters.molecular_search.score_method = "prob_score"
    MSParameters.molecular_search.output_score_method = "prob_score"

    print("\nLoading file: "+ file)

    parser = rawFileReader.ImportMassSpectraThermoMSFileReader(file)

    tic=parser.get_tic(ms_type='MS')[0]
    tic_df=pd.DataFrame({'time': tic.time,'scan': tic.scans})

    results = []
    for timestart in times:
        
        print('\nAveraging scans from %s-%s min' %(timestart, timestart+interval))
        scans=tic_df[tic_df.time.between(timestart,timestart+interval)].scan.tolist()
        
        mass_spectrum = parser.get_average_mass_spectrum_by_scanlist(scans)

        MzDomainCalibration(mass_spectrum, refmasslist,mzsegment=[0,1000]).run()

        mass_spectrum.molecular_search_settings.min_dbe = 0
        mass_spectrum.molecular_search_settings.max_dbe = 20

        mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1, 50)
        mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4, 100)
        mass_spectrum.molecular_search_settings.usedAtoms['O'] = (0, 15)
        mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0, 8)
        mass_spectrum.molecular_search_settings.usedAtoms['S'] = (0, 0)
        mass_spectrum.molecular_search_settings.usedAtoms['I'] = (0, 0)

        mass_spectrum.molecular_search_settings.isProtonated = True
        mass_spectrum.molecular_search_settings.isRadical = False
        mass_spectrum.molecular_search_settings.isAdduct = False

        mass_spectrum.molecular_search_settings.used_atom_valences = {'C': 4,
                                                                        '13C': 4,
                                                                        'H': 1,
                                                                        'O': 2,
                                                                        'N': 3,
                                                                        'S': 2,
                                                                        'I': 1
                                                                        }

        SearchMolecularFormulas(mass_spectrum, first_hit=True).run_worker_mass_spectrum()
        mass_spectrum.percentile_assigned(report_error=True)


        mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1, 50)
        mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4, 100)
        mass_spectrum.molecular_search_settings.usedAtoms['O'] = (0, 15)
        mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0, 8)
        mass_spectrum.molecular_search_settings.usedAtoms['S'] = (0, 1)
        mass_spectrum.molecular_search_settings.usedAtoms['I'] = (0, 1)

        mass_spectrum.molecular_search_settings.isProtonated = True
        mass_spectrum.molecular_search_settings.isRadical = False
        mass_spectrum.molecular_search_settings.isAdduct = False

        mass_spectrum.molecular_search_settings.used_atom_valences = {'C': 4,
                                                                        '13C': 4,
                                                                        'H': 1,
                                                                        'O': 2,
                                                                        'N': 3,
                                                                        'S': 2,
                                                                        'I': 1
                                                                        }
    
        SearchMolecularFormulas(mass_spectrum, first_hit=True).run_worker_mass_spectrum()
        mass_spectrum.percentile_assigned(report_error=True)



        
        assignments=mass_spectrum.to_dataframe()
        assignments['Time']=timestart
        results.append(assignments)
    
    results=pd.concat(results,ignore_index=True)
    
    return(results)



if __name__ == '__main__':

    data_dir = '/Users/christiandewey/data/Syn-cultures-May_June-2023/7803/'
    esifile = data_dir + 'Syn7803_neg.raw'
    icpmsfile = data_dir + 'Syn_7803.csv'
    refmasslist = '/Users/christiandewey/CoreMS/db/Hawkes_neg.ref'

    interval = 2
    '''#peak A
    time_min = 12
    time_max = 16

    #peak B
    time_min = 37
    time_max = 41'''

    #peak C
    time_min = 44
    time_max = 46

    '''#peak D
    time_min = 46
    time_max = 48

    #peak E
    time_min = 48
    time_max = 52'''


    os.chdir(data_dir)

    #icp_subset = subset_icpdata(icp_data_file = icpmsfile, heteroAtom = '127I', timerange = list(range(time_min,time_max,interval)), offset = -42)
    
    #icp_interp = interpolate(esi_file = esifile, icpsub=icp_subset, heteroAtom = '127I', timerange = list(range(time_min,time_max,interval)))

    assignment_results = assign_formula(file = esifile, times = list(range(time_min,time_max,interval)), interval=interval, refmasslist=refmasslist)

    assignment_results.to_csv(data_dir+'7803_peakC_neg_assignments.csv')

    #assignment_results = pd.read_csv(data_dir + '7803-neg_assignments.csv')

    #eic_dict, mass_spectrum = get_eics(esi_file=esifile, assignments = assignment_results, timerange = list(range(time_min,time_max,interval)))

    #results, bestresults = correlate(icp_interp=icp_interp, EICdic=eic_dict,heteroAtom='127I',assignments = assignment_results, timerange=list(range(time_min,time_max,interval)),threshold = 0.5)


    #os.system("sh /Users/christiandewey/CoreMS/reset_docker.sh")