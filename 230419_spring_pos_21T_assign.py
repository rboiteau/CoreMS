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
CoreMS run script for spring-env samples, collected at NHMFL in Nov 2023

Christian Dewey
18 April 2023
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

    print('\nAssignment took %.2f %s to complete' %(elapsed_time, unit))
    print('Started at ' + startdt_str )
    print('Finished at ' + enddt_str )


def assign_formula(esifile, times, charge, cal_ppm_threshold=(-1,1), refmasslist=None):
    
    corder=2

    parser = rawFileReader.ImportMassSpectraThermoMSFileReader(esifile)

    tic=parser.get_tic(ms_type='MS')[0]
    tic_df=pd.DataFrame({'time': tic.time,'scan': tic.scans})

    results=[]
    
    #calAssingmentParams()
    setAssingmentParams(charge)

    for timestart in times:
        print('\nfile: %s\ntimestart:%s'  %(esifile,timestart))
        scans=tic_df[tic_df.time.between(timestart,timestart+interval)].scan.tolist()

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

            pmzrfs = pd.DataFrame(mzrefs)
            pmzrfs.to_csv('cal_mzs_%s.csv' %esifile.split('.')[0])
            calfn.recalibrate_mass_spectrum(mass_spectrum, imzmeas, mzrefs, order=corder)

        print('\nassigning with first parameter set...')

        SearchMolecularFormulas(mass_spectrum, first_hit=False).run_worker_mass_spectrum()

        #print('\nresults with first parameter set...')
        mass_spectrum.percentile_assigned(report_error=True)

        #print('\nassigning with second parameter set...')
        #setAssingmentParams2(charge)

        #SearchMolecularFormulas(mass_spectrum, first_hit=True).run_worker_mass_spectrum()

        #print('\nresults with second parameter set...')
        #mass_spectrum.percentile_assigned(report_error=True)

        assignments=mass_spectrum.to_dataframe()

        assignments['Time']=timestart

        results.append(assignments)
    
    results=pd.concat(results,ignore_index=True)

    return(results)   


def calAssingmentParams():
    # set assignment parameters
    MSParameters.mass_spectrum.threshold_method = 'signal_noise'
    MSParameters.mass_spectrum.s2n_threshold = 3

    MSParameters.molecular_search.error_method = 'None'
    MSParameters.molecular_search.min_ppm_error = -0.25
    MSParameters.molecular_search.max_ppm_error = 0.25

    MSParameters.molecular_search.isProtonated = True
    MSParameters.molecular_search.isRadical = False
    MSParameters.molecular_search.isAdduct = False

    MSParameters.molecular_search.score_method = "prob_score"
    MSParameters.molecular_search.output_score_method = "prob_score"

    MSParameters.molecular_search.url_database = 'postgresql+psycopg2://coremsappdb:coremsapppnnl@localhost:5432/coremsapp'
    MSParameters.molecular_search.min_dbe = -1
    MSParameters.molecular_search.max_dbe = 20
    MSParameters.molecular_search.ion_charge = 1 # absolute value; multiplied by polarity w/in code

    MSParameters.molecular_search.usedAtoms['C'] = (1,50)  
    MSParameters.molecular_search.usedAtoms['H'] = (4,100)
    MSParameters.molecular_search.usedAtoms['O'] = (0,10)
    MSParameters.molecular_search.usedAtoms['N'] = (0,3)
    MSParameters.molecular_search.usedAtoms['S'] = (0,0)
    MSParameters.molecular_search.usedAtoms['P'] = (0,0)
    MSParameters.molecular_search.usedAtoms['Na'] = (0,0)
    MSParameters.molecular_search.usedAtoms['Cu'] = (0,0)
    MSParameters.molecular_search.usedAtoms['K'] = (0,0)
    MSParameters.molecular_search.usedAtoms['Fe'] = (0,0)


def setAssingmentParams(ion_charge):
    # set assignment parameters
    MSParameters.mass_spectrum.threshold_method = 'signal_noise'
    MSParameters.mass_spectrum.s2n_threshold = 3

    MSParameters.molecular_search.error_method = 'None'
    MSParameters.molecular_search.min_ppm_error = -0.25
    MSParameters.molecular_search.max_ppm_error = 0.25

    MSParameters.molecular_search.isProtonated = True
    MSParameters.molecular_search.isRadical = False
    MSParameters.molecular_search.isAdduct = False

    MSParameters.molecular_search.score_method = "prob_score"
    MSParameters.molecular_search.output_score_method = "prob_score"

    MSParameters.molecular_search.url_database = 'postgresql+psycopg2://coremsappdb:coremsapppnnl@localhost:5432/coremsapp'
    MSParameters.molecular_search.min_dbe = -1
    MSParameters.molecular_search.max_dbe = 20
    MSParameters.molecular_search.ion_charge = ion_charge # absolute value; multiplied by polarity w/in code

    MSParameters.molecular_search.usedAtoms['C'] = (1,100)  
    MSParameters.molecular_search.usedAtoms['H'] = (4,200)
    MSParameters.molecular_search.usedAtoms['O'] = (0,20)
    MSParameters.molecular_search.usedAtoms['N'] = (0,20)
    MSParameters.molecular_search.usedAtoms['S'] = (0,4)
    MSParameters.molecular_search.usedAtoms['P'] = (0,4)
    MSParameters.molecular_search.usedAtoms['Na'] = (0,1)
    MSParameters.molecular_search.usedAtoms['Cu'] = (0,0)
    MSParameters.molecular_search.usedAtoms['K'] = (0,0)
    MSParameters.molecular_search.usedAtoms['Fe'] = (0,0)


def setAssingmentParams2(ion_charge):
    # set assignment parameters
    MSParameters.mass_spectrum.threshold_method = 'signal_noise'
    MSParameters.mass_spectrum.s2n_threshold = 3

    MSParameters.molecular_search.error_method = 'None'
    MSParameters.molecular_search.min_ppm_error = -0.25
    MSParameters.molecular_search.max_ppm_error = 0.25

    MSParameters.molecular_search.isProtonated = True
    MSParameters.molecular_search.isRadical = False
    MSParameters.molecular_search.isAdduct = False

    MSParameters.molecular_search.score_method = "prob_score"
    MSParameters.molecular_search.output_score_method = "prob_score"

    MSParameters.molecular_search.url_database = 'postgresql+psycopg2://coremsappdb:coremsapppnnl@localhost:5432/coremsapp'
    MSParameters.molecular_search.min_dbe = -1
    MSParameters.molecular_search.max_dbe = 20
    MSParameters.molecular_search.ion_charge = ion_charge # absolute value; multiplied by polarity w/in code

    MSParameters.molecular_search.usedAtoms['C'] = (10,100)  
    MSParameters.molecular_search.usedAtoms['H'] = (20,150)
    MSParameters.molecular_search.usedAtoms['O'] = (0,20)
    MSParameters.molecular_search.usedAtoms['N'] = (0,20)
    MSParameters.molecular_search.usedAtoms['S'] = (0,2)
    MSParameters.molecular_search.usedAtoms['P'] = (0,2)
    MSParameters.molecular_search.usedAtoms['Na'] = (0,2)
    MSParameters.molecular_search.usedAtoms['Cu'] = (0,1)
    MSParameters.molecular_search.usedAtoms['K'] = (0,1)
    MSParameters.molecular_search.usedAtoms['Fe'] = (0,1)


if __name__ == '__main__':

    start = time.time()  #for duration
    startdt = datetime.now()
    
    data_dir = '/mnt/disks/orca-data/mz-windowing/pos/spring/'

    fname = '230419_spring-env_pos_001.csv'

    mzref = "/home/CoreMS/tests/tests_data/ftms/nom_pos.ref" 
    
    interval = 4     # window in which scans are averaged
    time_range = [8,12]    

    results = []
    
    times = list(range(time_range[0],time_range[1],interval))

    flist=os.listdir(data_dir)
    f_raw = [f for f in flist if '.raw' in f]
    os.chdir(data_dir)
    
    ion_charge = 1
    i = 1
    for f in f_raw:
        if 'spring_fullmz_rep2' in f:
            print("\n\n\n%s/%s files" %(i, len(f_raw)))
            output = assign_formula(esifile = f, times = times, charge = ion_charge, cal_ppm_threshold=(-1,1), refmasslist = mzref)
            output['file'] = f 
            results.append(output)
            i = i + 1
    
    df = pd.concat(results)
    df.to_csv(data_dir+fname)

    printRunTime()

    
