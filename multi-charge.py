# Import the os module
import os
import warnings
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')
import sys
sys.path.append('./')

from corems.mass_spectra.input import rawFileReader
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.encapsulation.factory.parameters import MSParameters
from corems.mass_spectrum.calc.Calibration import MzDomainCalibration


def assign_formula(file, times, interval): 
    MSParameters.molecular_search.error_method = 'None'
    MSParameters.molecular_search.min_ppm_error = -0.25
    MSParameters.molecular_search.max_ppm_error = 0.25
    MSParameters.molecular_search.ion_charge = 2
    MSParameters.molecular_search.db_chunk_size = 500

    MSParameters.mass_spectrum.min_calib_ppm_error = -1
    MSParameters.mass_spectrum.max_calib_ppm_error = 1
    MSParameters.mass_spectrum.calib_pol_order = 2
    #MSParameters.mass_spectrum.calib_sn_threshold = 3
    MSParameters.mass_spectrum.min_picking_mz=100
    MSParameters.mass_spectrum.max_picking_mz=900
    MSParameters.mass_spectrum.threshold_method = 'log'
    MSParameters.mass_spectrum.log_nsigma=400

    MSParameters.ms_peak.peak_min_prominence_percent = 0.001

    MSParameters.molecular_search.url_database = 'postgresql+psycopg2://coremsappdb:coremsapppnnl@localhost:5432/coremsapp' #'postgresql+psycopg2://coremsappdb:coremsapppnnl@molformdb-1:5432/coremsapp'
    MSParameters.molecular_search.score_method = "prob_score"
    MSParameters.molecular_search.output_score_method = "prob_score"


    print("Loading file: "+ file)
    # Read in sample list and load MS data
    MSfiles={}
    parser = rawFileReader.ImportMassSpectraThermoMSFileReader(file)

    MSfiles[file]=parser

    tic=parser.get_tic(ms_type='MS')[0]
    tic_df=pd.DataFrame({'time': tic.time,'scan': tic.scans})

    refmasslist = "/Users/christiandewey/CoreMS/db/nom_pos.ref"

    results = []
    for timestart in times:
        print(timestart)
        scans=tic_df[tic_df.time.between(timestart,timestart+interval)].scan.tolist()
        
        mass_spectrum = parser.get_average_mass_spectrum_by_scanlist(scans)

        MzDomainCalibration(mass_spectrum, refmasslist,mzsegment=[0,1000]).run()

        #First assignment iteration (CHON with adducts)
        mass_spectrum.molecular_search_settings.min_dbe = 0
        mass_spectrum.molecular_search_settings.max_dbe = 30

        mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1, 70)
        mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4, 100)
        mass_spectrum.molecular_search_settings.usedAtoms['O'] = (0, 20)
        mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0, 10)
        mass_spectrum.molecular_search_settings.usedAtoms['P'] = (0, 1)
        mass_spectrum.molecular_search_settings.usedAtoms['Fe'] = (0, 1)
        mass_spectrum.molecular_search_settings.usedAtoms['Cl'] = (0, 1)

        mass_spectrum.molecular_search_settings.isProtonated = True
        mass_spectrum.molecular_search_settings.isRadical = False
        mass_spectrum.molecular_search_settings.isAdduct = False
        #mass_spectrum.molecular_search_settings.max_oc_filter=1.2
        #mass_spectrum.molecular_search_settings.max_hc_filter=3
        mass_spectrum.molecular_search_settings.used_atom_valences = {'C': 4,
                                                                        '13C': 4,
                                                                        'H': 1,
                                                                        'O': 2,
                                                                        'N': 3,
                                                                        'P': 3,
                                                                        'Fe': 3,
                                                                        'Cl': 1
                                                                        }

    
        SearchMolecularFormulas(mass_spectrum, first_hit=True).run_worker_mass_spectrum()
        mass_spectrum.percentile_assigned(report_error=True)
        
        assignments=mass_spectrum.to_dataframe()
        assignments['Time']=timestart
        results.append(assignments)
    
    results=pd.concat(results,ignore_index=True)
    
    return(results)



if __name__ == '__main__':

    data_dir = '/Users/christiandewey/data-temp/'
    results = []

    interval = 2
    time_min = 10
    time_max = 16
    times = list(range(time_min,time_max,interval))

    flist = os.listdir(data_dir)
    f_raw = [f for f in flist if '.raw' in f]
    os.chdir(data_dir)
    i=1

    for f in f_raw:
        print(f)
        output = assign_formula(file = f, times = times, interval=interval)
        output['file'] = f
        #make_plots(output,f)
        fname = 'test_assignments.csv'
        output.to_csv(data_dir+fname)
        i = i + 1 

    #os.system("sh /Users/christiandewey/CoreMS/reset_docker.sh")
