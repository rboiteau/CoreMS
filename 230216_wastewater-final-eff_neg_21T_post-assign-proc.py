import os
from tempfile import tempdir
import time
import numpy as np
import warnings
from datetime import date, datetime
import pandas as pd

warnings.filterwarnings("ignore")
from pathlib import Path
import sys
sys.path.append("./")

from corems.mass_spectra.input import rawFileReader
from corems.molecular_id.factory.classification import HeteroatomsClassification, Labels
from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.encapsulation.factory.parameters import MSParameters
from corems.encapsulation.constant import Atoms
from corems.mass_spectrum.calc.Calibration import MzDomainCalibration
import corems.lc_icpms_ftms.calc.lc_icrms_qc_assign as icrms
import corems.lc_icpms_ftms.calc.lc_icrms_helpers as lcmsfns

"""
Post-assignment processing for m/z windowing project

Intended for wastewater-final-eff_neg_21TNHMFL-Jan_001.csv


Adds following data columns:

    mol_class, Window Size (m/z), m/z window, Rep


Performs blank subtraction

Christian Dewey
16 Feb 23
"""

def printRunTime():

    global start, startdt

    elapsed_time_sec = (time.time() - start) 

    if elapsed_time_sec > 3600:
        elapsed_time = elapsed_time_sec / 3600
        unit = 'hr'
    else:
        elapsed_time = elapsed_time_sec / 60
        unit = 'min'

    enddt = datetime.now()

    startdt_str = startdt.strftime("%d-%b-%Y %H:%M:%S %z")
    enddt_str = enddt.strftime("%d-%b-%Y %H:%M:%S %z")

    print('\nRan in %.2f %s' %(elapsed_time, unit))
    print('Started at ' + startdt_str )
    print('Finished at ' + enddt_str )



def postAssignProcess(heter):

    ##add 'Window Size (m/z)', 'm/z window', 'Rep', 'mol_class' columns

    global results

    results = lcmsfns.add_mzwindow_col(results)                 # adds 'Window Size (m/z)' and 'm/z window' columns
    results = lcmsfns.addRepCol(results)                        # adds 'Rep' column
    molclasses = lcmsfns.get_mol_class(heter)                   # creates list of mol_classes based on heteroatom list
    results = lcmsfns.assign_mol_class(results,molclasses)      # adds 'mol_class' column 



    ## adds column with blank file identity to use in blank subtraction  

    raw_filelist = results['file'].unique()
    blank_files = [f for f in raw_filelist if 'qh2o' in f]
    blank_data = []

    for f in blank_files:
            
        temp = results[results['file'] == f] 

        blank_data.append(temp)

    blanks_df = pd.concat(blank_data)  

    results['blank file'] = results.index

    for window in results['m/z window'].unique():

        temp1 = results[results['m/z window'] == window]        # all features collected in given m/z window       
        btemp1 = blanks_df[blanks_df['m/z window'] == window]   # all blank feautres collected in same m/z window
        
        for r in temp1['Rep'].unique():

            temp2 = temp1[temp1['Rep'] == r]
            btemp2 = btemp1[btemp1['Rep'] == r]

            temp2['blank file'] = btemp2['file'].iloc[0]
            results[(results['m/z window'] == window) & (results['Rep'] == r)]  = temp2





if __name__ == '__main__':
    start = time.time()  #for duration
    startdt = datetime.now()

    asgn_dir = '/Users/christiandewey/Library/CloudStorage/GoogleDrive-christian.w.dewey@gmail.com/My Drive/manuscripts/2023_Dewey-Boiteau-etal_mz-windowing/assignments/wastewater-final-eff/'
    fname = '230215_wastewater-final-eff_neg_21TNHMFL-Jan_001.csv'

    heter = ['N',  'S', 'P']  # from 230215_wastewater-final-eff_neg_21T_assign.py


    results = pd.read_csv(asgn_dir+fname)



    postAssignProcess(heter)

    results.to_csv(asgn_dir+'p'+fname)

    printRunTime()

    
