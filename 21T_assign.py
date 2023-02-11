import os
from tempfile import tempdir
import time
from turtle import color
import numpy as np
import warnings
from datetime import date, datetime

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
import pandas as pd


def setAssingmentParams():
    # set assignment parameters
    MSParameters.mass_spectrum.threshold_method = 'signal_noise'
    MSParameters.mass_spectrum.s2n_threshold = 2
    MSParameters.ms_peak.peak_min_prominence_percent = 0.001

    MSParameters.molecular_search.error_method = 'None'
    MSParameters.molecular_search.min_ppm_error = -0.25
    MSParameters.molecular_search.max_ppm_error = 0.25

    MSParameters.molecular_search.isProtonated = True
    MSParameters.molecular_search.isRadical = False
    MSParameters.molecular_search.isAdduct = False

    MSParameters.molecular_search.score_method = "prob_score"
    MSParameters.molecular_search.output_score_method = "prob_score"

    MSParameters.molecular_search.url_database = None
    MSParameters.molecular_search.min_dbe = -1
    MSParameters.molecular_search.max_dbe = 20

    MSParameters.molecular_search.usedAtoms['C'] = (1,40)
    MSParameters.molecular_search.usedAtoms['H'] = (4,80)
    MSParameters.molecular_search.usedAtoms['O'] = (1,18)
    MSParameters.molecular_search.usedAtoms['N'] = (0,4)
    MSParameters.molecular_search.usedAtoms['Fe'] = (0,1)
    MSParameters.molecular_search.usedAtoms['S'] = (0,2)
    MSParameters.molecular_search.usedAtoms['P'] = (0,2)
    MSParameters.molecular_search.usedAtoms['Na'] = (0,1)
    MSParameters.molecular_search.usedAtoms['Cu'] = (0,1)

def assign_formula(esifile, times, refmasslist=None):

    cal_ppm_threshold=(-1,1)
    charge=1
    corder=2

    parser = rawFileReader.ImportMassSpectraThermoMSFileReader(esifile)

    tic=parser.get_tic(ms_type='MS')[0]
    tic_df=pd.DataFrame({'time': tic.time,'scan': tic.scans})

    results=[]
    
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

            calfn.recalibrate_mass_spectrum(mass_spectrum, imzmeas, mzrefs, order=corder)


        SearchMolecularFormulas(mass_spectrum, first_hit=False).run_worker_mass_spectrum()

        mass_spectrum.percentile_assigned(report_error=True)

        assignments=mass_spectrum.to_dataframe()

        assignments['Time']=timestart

        results.append(assignments)
    
    results=pd.concat(results,ignore_index=True)

    return(results)   





if __name__ == '__main__':
    
    data_dir = '/Volumes/Samsung_T5/NHMFL/clustering/pos/spring/'
    svdir = '/Users/christiandewey/Library/CloudStorage/GoogleDrive-christian.w.dewey@gmail.com/My Drive/manuscripts/2023_Dewey-Boiteau-etal_mz-windowing/assignments/'
    fname = 'spring_21TNHMFL_nov22-run_001.csv'
    mzref = "/Users/christiandewey/CoreMS/tests/tests_data/ftms/nom_pos.ref"

    interval = 2
    time_range = [4,24]
    times = list(range(time_range[0],time_range[1],interval))

    flist=os.listdir(data_dir)
    os.chdir(data_dir)

    setAssingmentParams()
    results = []
    start = time.time()  #for loop duration

    startdt = datetime.now()
    

    f_raw = [f for f in flist if '.raw' in f]
    i = 1
    for f in f_raw:
        print("\n\n\n%s/%s files" %(i, len(f_raw)))
        output = assign_formula(f,times,mzref)
        output['file'] = f 
        results.append(output)
        i = i + 1

    df = pd.concat(results)
    df.to_csv(svdir+fname)

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