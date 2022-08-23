#### Matches LC-ICPMS data and LC-ESIMS data and assignes MF to peak.

#### Christian Dewey
#### 18 August 2022

import os
import pandas as pd
import numpy as np
import warnings
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

class Aligned_ICP_ESI:
    """
    classdocs
    """

    def __init__(self, icp_file_location,esi_file_location):
        self._icpfile = icp_file_location
        self._esifile = esi_file_location
        self.esi_parser = None 
        self.icp_parser = None 
        self.offset = -27
        self.timestart = 8.2 # min
        self.timestop = 9.3 # min
        #self.elements = ['63Cu']
        #self.etime = ['Time ' + m for m in self.elements]
        self.element = '63Cu'
        self.etime = 'Time ' + self.element
        self.threshold = 0.8

    def getMSData(self):    
        MSParameters.mass_spectrum.threshold_method = 'signal_noise'
        MSParameters.mass_spectrum.s2n_threshold = 6
        MSParameters.ms_peak.peak_min_prominence_percent = 0.1
        self.esi_parser = rawFileReader.ImportMassSpectraThermoMSFileReader(self._esifile)
        self.icp_data = pd.read_csv(self._icpfile)

    def subset_icpdata(self):
        element = self.element
        etime = self.etime
        icp_subset = self.icp_data[[element,etime]]
        icp_subset[etime] = (icp_subset[etime] + self.offset) / 60 
        icp_subset = icp_subset[icp_subset[etime].between(self.timestart,self.timestop)]
        icp_subset.plot.line(x=etime,y=element)
        #plt.show()
        return icp_subset
    
    def subset_esidata(self):
        element = self.element
        etime = self.etime
        tic = self.esi_parser.get_tic(ms_type = 'MS')[0]
        tic_df=pd.DataFrame({'time': tic.time,'scan': tic.scans})
        scans=tic_df[tic_df.time.between(self.timestart,self.timestop)].scan.tolist()
        AverageMS = self.esi_parser.get_average_mass_spectrum_by_scanlist(scans)
        #AverageMS.plot_mz_domain_profile()
        #plt.show()
        ### print(AverageMS.mz_exp.size)
        ### The correct assignment should be 381m/z

        EICdic = {}
        pbar = tqdm.tqdm(AverageMS.mz_exp, desc="Getting EICs")
        for mz in pbar:
            EIC=self.esi_parser.get_eics(target_mzs=[mz],tic_data={},peak_detection=False,smooth=False)
            EICdic[mz]=EIC[0][mz]

        ###Interpolate LC-ICPMS data to obtain times matching ESI data
        times=tic_df[tic_df.time.between(self.timestart,self.timestop)].time.tolist()
        icpsubset2 = self.subset_icpdata()
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

        mscorr={}
        EICcurr=pd.DataFrame(index=icp_interp[etime],columns=['EIC',element])

        EICcurr=pd.DataFrame(index=icp_interp[etime],columns=['EIC',element])
        EICcurr[element]=icp_interp[element].array

        pbar = tqdm.tqdm(EICdic.keys(),desc="Running correlation")
        for mz in pbar:
            EIC=pd.DataFrame({'EIC':EICdic[mz].eic,'Time':EICdic[mz].time})
            EIC_sub=EIC[EIC['Time'].between(self.timestart,self.timestop)]
            #print(EIC_sub['EIC'])

            EICcurr['EIC']=EIC_sub['EIC'].values

            corvalue=EICcurr.corr(method='pearson')
            mscorr[mz]=corvalue.EIC[element]**2

        mzs_corr = pd.DataFrame.from_dict(mscorr,orient='index',columns=['corr'])
        
        return mzs_corr, AverageMS

    def assignFormulas(self):
        threshold = self.threshold
        mzs_corr, mass_spectrum = self.subset_esidata()

        #Get molecular formula of average mass spectrum. 

        mass_spectrum.molecular_search_settings.error_method = 'None'
        mass_spectrum.molecular_search_settings.min_ppm_error = -2
        mass_spectrum.molecular_search_settings.max_ppm_error = 2

        mass_spectrum.molecular_search_settings.url_database = None
        mass_spectrum.molecular_search_settings.min_dbe = 0
        mass_spectrum.molecular_search_settings.max_dbe = 20

        mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1, 50)
        mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4, 100)
        mass_spectrum.molecular_search_settings.usedAtoms['O'] = (1, 20)
        mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0, 4)
        mass_spectrum.molecular_search_settings.usedAtoms['S'] = (0, 0)
        mass_spectrum.molecular_search_settings.usedAtoms['Cl'] = (0, 0)
        mass_spectrum.molecular_search_settings.usedAtoms['Br'] = (0, 0)
        mass_spectrum.molecular_search_settings.usedAtoms['P'] = (0, 0)
        mass_spectrum.molecular_search_settings.usedAtoms['Na'] = (0, 0)
        mass_spectrum.molecular_search_settings.usedAtoms['Cu'] = (0,1)

        mass_spectrum.molecular_search_settings.isProtonated = True
        mass_spectrum.molecular_search_settings.isRadical = False
        mass_spectrum.molecular_search_settings.isAdduct = False

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

        assignments.insert(4,'corr',mzs_corr['corr'].values)

        threshold=0.8
        results=assignments[assignments['corr']>threshold].filter(['m/z','corr','Peak Height','Confidence Score','Molecular Formula','Cu'])
        bestresults=results[results['Cu']==1]
        return bestresults