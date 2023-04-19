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

Intended for 230301_spring-env_pos_21TNHMFL-Nov-22_002.csv


Adds following data columns:

    mol_class, Window Size (m/z), m/z window, Rep, blank file

Christian Dewey
01 Mar 23
"""



def postAssignProcess(heter):

    ##add 'Window Size (m/z)', 'm/z window', 'Rep', 'mol_class' columns

    global results

    results = lcmsfns.add_mzwindow_col(results)                 # adds 'Window Size (m/z)' and 'm/z window' columns
    results = lcmsfns.addRepCol(results)                        # adds 'Rep' column
    molclasses = lcmsfns.get_mol_class(heter)                   # creates list of mol_classes based on heteroatom list
    results = lcmsfns.assign_mol_class(results,molclasses)      # adds 'mol_class' column 



def blankFileCreate():

    # create 200 m/z blank files from 100 m/z blank files and add column with blank file identity 
    global results 
    
    flist = results.file.unique()
    blank_files = [f for f in flist if 'qh2o' in f]

    blank_data = []

    for f in blank_files:

        if 'fullmz' not in f:
            
            temp = results[results['file'] == f] 

            blank_data.append(temp)

    blanks_df = pd.concat(blank_data)

    rep1_temp = blanks_df[~blanks_df['file'].str.contains('rep2')]
    rep1_temp = rep1_temp[rep1_temp['m/z'] <= 600]
    rep1_temp['file'] = 'mz200_400_600_blnk'
    rep1_temp['m/z window'] = '400-600 m/z'

    blanks_df = pd.concat([blanks_df,rep1_temp])

    rep2_temp = blanks_df[blanks_df['file'].str.contains('rep2')]
    rep2_temp = rep2_temp[rep2_temp['m/z'] <= 600]
    rep2_temp['file'] = 'mz200_400_600_blnk_rep2'
    rep2_temp['m/z window'] = '400-600 m/z'


    blanks_df = pd.concat([blanks_df,rep2_temp])

    print(blanks_df['file'].unique())

    print(blanks_df['m/z window'].unique())

    mz200_blanks = blanks_df[blanks_df['m/z window'] == '400-600 m/z']
    mz200_blanks['Window Size (m/z)'] = '200'

    n200bs = len(mz200_blanks)
    max_index = max(results.index)
    compatible_index = range(max_index, max_index+n200bs,1)
    mz200_blanks.index = compatible_index

    results = pd.concat([results, mz200_blanks])

    blank_data = []

    for f in blank_files:

        if 'fullmz' in f:
            
            temp = results[results['file'] == f] 

            blank_data.append(temp)

    blanks_fullmz_df = pd.concat(blank_data)
    blanks_df = pd.concat([blanks_df,blanks_fullmz_df])

    df_bs = []
    for window in results['m/z window'].unique():

        temp1 = results[results['m/z window'] == window] # all features collected in given m/z window       
        btemp1 = blanks_df[blanks_df['m/z window'] == window] # all blank feautres collected in same m/z window
        for r in temp1['Rep'].unique():
            temp2 = temp1[temp1['Rep'] == r]
            btemp2 = btemp1[btemp1['Rep'] == r]

            temp2['blank file'] = btemp2['file'].iloc[0]

            df_bs.append(temp2)

    df_bs = pd.concat(df_bs)

    results =df_bs



if __name__ == '__main__':

    data_dir = '/mnt/disks/orca-data/mz-windowing/pos/spring/'

    flist = [f for f in os.listdir() if '.csv' in f]

    for fname in flist:

        results = pd.read_csv(data_dir+fname)

        heter = ['N', 'Na', 'S', 'P', 'K', 'Cu','Fe']  # from 2300301_spring-env_pos_21T_assign.py

        postAssignProcess(heter)

        blankFileCreate()

        results.to_csv(data_dir+'p'+fname)

    
