#### ICPMS ESIMS alignment helpers

#### Christian Dewey
#### 23 January 2023

import os
import pandas as pd
import numpy as np
import warnings
import re
warnings.filterwarnings("ignore")

import sys
sys.path.append("./")

from pathlib import Path

import tqdm

import matplotlib.pyplot as plt
# from PySide2.QtWidgets import QFileDialog, QApplication
# from PySide2.QtCore import Qt

from corems.mass_spectra.input import rawFileReader
from corems.molecular_id.factory.classification import HeteroatomsClassification, Labels
from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.encapsulation.factory.parameters import MSParameters

            

def subset_icpdata(icp_data = None, heteroAtom = '127I', timerange = [0,1], offset = 0):
    timestart = timerange[0]
    timestop = timerange[1]
    element = heteroAtom
    etime = 'Time ' + heteroAtom
    icp_subset = icp_data[[element,etime]]
    icp_subset[etime] = (icp_subset[etime] + offset) / 60 
    icp_subset = icp_subset[icp_subset[etime].between(timestart,timestop)]

    return icp_subset

def get_eics(esi_parser=None, timerange = [0,1]):

    timestart = timerange[0]
    timestop = timerange[1]
    tic = esi_parser.get_tic(ms_type = 'MS')[0]
    tic_df=pd.DataFrame({'time': tic.time,'scan': tic.scans})
    scans=tic_df[tic_df.time.between(timestart,timestop)].scan.tolist()
    AverageMS = esi_parser.get_average_mass_spectrum_by_scanlist(scans)

    EICdic = {}
    pbar = tqdm.tqdm(AverageMS.mz_exp, desc="Getting EICs")
    
    for mz in pbar:   
        #   print('AverageMS mz:' + str(mz))
        EIC=esi_parser.get_eics(target_mzs=[mz],tic_data={},peak_detection=False,smooth=False)
        EICdic[mz]=EIC[0][mz]

    return EICdic, AverageMS

def interp(esi_parser = None, heteroAtom = '127I', timerange = [0,1], icpsub=None):
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
        holder[row,1] = mzs_corr[mzs_corr.index == mz].iloc[0]['corr']
        holder[row,0] = mz

    pdholder = pd.DataFrame(holder, columns = ['m/z', 'corr'])

    assignments.insert(4,'corr',pdholder['corr'].values)

    assignments.insert(5,'m/z-2',pdholder['m/z'].values)

    match = re.findall(r'[A-Za-z]+|\d+', heteroAtom)
    heteroAtom = match[1]
    
    results=assignments[assignments['corr']>threshold].filter(['m/z','corr','Peak Height','Confidence Score','Molecular Formula',heteroAtom])

    bestresults=results[(results[heteroAtom]>=1)]

    return results, bestresults

"""
    def assignFormulas(self,  threshold):
        # elementDict = {'C': (1,50), 'H':(4,100), etc}
       # threshold = self.threshold
        print('threshold: ' + str(threshold))
        print('offset: ' + str(self.offset))
        mzs_corr, mass_spectrum = self.subset_esidata()

        #Get molecular formula of average mass spectrum. 

        '''        mass_spectrum.molecular_search_settings.error_method = 'None'
        mass_spectrum.molecular_search_settings.min_ppm_error = -.25
        mass_spectrum.molecular_search_settings.max_ppm_error = 0.25

        mass_spectrum.molecular_search_settings.url_database = None
        mass_spectrum.molecular_search_settings.min_dbe = 0
        mass_spectrum.molecular_search_settings.max_dbe = 100

        self.elementDict = elementDict
        elementList = elementDict.keys()

        for e in elementList:
            mass_spectrum.molecular_search_settings.usedAtoms[e] = self.elementDict[e]


        mass_spectrum.molecular_search_settings.isProtonated = True
        mass_spectrum.molecular_search_settings.isRadical = False
        mass_spectrum.molecular_search_settings.isAdduct = False'''

        # mass_spectrum.filter_by_max_resolving_power(15, 2)
        SearchMolecularFormulas(mass_spectrum, first_hit=False).run_worker_mass_spectrum()

        mass_spectrum.percentile_assigned(report_error=True)
        mass_spectrum.molecular_search_settings.score_method = "prob_score"
        mass_spectrum.molecular_search_settings.output_score_method = "prob_score"

        mass_spectrum_by_classes = HeteroatomsClassification(mass_spectrum, choose_molecular_formula=True)
        try:
            mass_spectrum_by_classes.plot_ms_assigned_unassigned()
            plt.show()
        except (RuntimeError, TypeError, NameError):
            pass
        
        try:
            mass_spectrum_by_classes.plot_mz_error()
            plt.show()
        except (RuntimeError, TypeError, NameError):
            pass       
    
        try:
            mass_spectrum_by_classes.plot_ms_class("O2")
            plt.show()
        except (RuntimeError, TypeError, NameError):
            pass

        assignments=mass_spectrum.to_dataframe()
        assignments=assignments.sort_values(by=['m/z'])


        holder = pd.DataFrame( index = range(len(assignments['m/z'])), columns = ['m/z', 'corr'])  #index = assignments['Index'],

        holder = np.zeros((len(assignments['m/z']),2))


        for mz,row in zip(assignments['m/z'], range(len(assignments['m/z']))):

            holder[row,1] = mzs_corr[mzs_corr.index == mz].iloc[0]['corr']
            holder[row,0] = mz
 
        pdholder = pd.DataFrame(holder, columns = ['m/z', 'corr'])

       # for i, j, mzi, mzj in zip(assignments['Index'],pdholder.index,assignments['m/z'], pdholder['m/z']):
        #    print(i, j, mzi, mzj)
        pdholder.to_csv('/Users/christiandewey/Downloads/mzs_corr.csv')

        

       # print('shape mzs_corr: ')
       # print(np.shape(holder))
       # print(holder)


        #print('shape assignments: ')
        ##print(np.shape(assignments))
        #print(assignments)
        


        assignments.insert(4,'corr',pdholder['corr'].values)

        assignments.insert(5,'mz2',pdholder['m/z'].values)

        assignments.to_csv('/Users/christiandewey/Downloads/assignments.csv')

       
        match = re.findall(r'[A-Za-z]+|\d+', self.heteroAtom)
        heteroAtom = match[1]
        
        
        results=assignments[assignments['corr']>threshold].filter(['m/z','corr','Peak Height','Confidence Score','Molecular Formula',heteroAtom])

        print(results)

        bestresults=results[(results[heteroAtom]>=1)]


        #print(bestresults)
        return bestresults


"""