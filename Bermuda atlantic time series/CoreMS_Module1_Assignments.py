# LC-MS formula assignment module
# This file should be run from the CoreMS main folder. 
# Data files and the sample lists should be placed in the CoreMS/usrdata subfolder. 


# Developed for data collected on DOM isolated from the Bermuda Atlantic Time Series
# RMB update 7/10/2023
# Contributors: Yuri Corilo, Will Kew, Christian Dewey, Rene Boiteau


# Import the os module
import os
from tempfile import tempdir
import warnings
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')
import sys
sys.path.append('./')
#os.chdir('/CoreMS')
from corems.mass_spectra.input import rawFileReader
from corems.molecular_id.factory.classification import HeteroatomsClassification, Labels
from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.encapsulation.factory.parameters import MSParameters
from corems.encapsulation.constant import Atoms
from corems.mass_spectrum.calc.Calibration import MzDomainCalibration
#os.chdir('/CoreMS/usrdata')
os.chdir('C://Users/boiteaur/Desktop/corems/rawfiles_bats')
#Function that averages and calibrates mass spectra and assigns molecular formula

def assign_formula(file, times):

    MSParameters.molecular_search.error_method = 'None'
    MSParameters.molecular_search.min_ppm_error = -0.3
    MSParameters.molecular_search.max_ppm_error = 0.3
    MSParameters.molecular_search.default_ion_charge = 1
    MSParameters.molecular_search.min_ion_charge = 1
    MSParameters.molecular_search.max_ion_charge = 2
    MSParameters.molecular_search.db_chunk_size = 500
    MSParameters.mass_spectrum.min_calib_ppm_error = 0.5
    MSParameters.mass_spectrum.max_calib_ppm_error = 2
    MSParameters.mass_spectrum.calib_pol_order = 2
    #MSParameters.mass_spectrum.calib_sn_threshold = 3
    MSParameters.mass_spectrum.min_picking_mz=220.001
    MSParameters.mass_spectrum.max_picking_mz=900
    MSParameters.mass_spectrum.threshold_method = 'signal_noise'
    MSParameters.mass_spectrum.s2n_threshold=2

    MSParameters.ms_peak.peak_min_prominence_percent = 0.01

    MSParameters.molecular_search.url_database = 'postgresql+psycopg2://coremsappdb:coremsapppnnl@corems-molformdb-1:5432/coremsapp'
    MSParameters.molecular_search.score_method = "prob_score"
    MSParameters.molecular_search.output_score_method = "prob_score"

    refmasslist = "Seawater_NOM_pos.ref"


    print("Loading file: "+ file)
    # Read in sample list and load MS data
    MSfiles={}
    parser = rawFileReader.ImportMassSpectraThermoMSFileReader(file)

    MSfiles[file]=parser

    tic=parser.get_tic(ms_type='MS')[0]
    tic_df=pd.DataFrame({'time': tic.time,'scan': tic.scans})


    results = []
    for timestart in times:
        print(timestart)
        scans=tic_df[tic_df.time.between(timestart,timestart+interval)].scan.tolist()

        mass_spectrum = parser.get_average_mass_spectrum_by_scanlist(scans)

        MzDomainCalibration(mass_spectrum, refmasslist,mzsegment=[0,1000]).run()

        #First assignment iteration (CHON with adducts)
        mass_spectrum.molecular_search_settings.min_dbe = 0
        mass_spectrum.molecular_search_settings.max_dbe = 20

        mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1, 50)
        mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4, 100)
        mass_spectrum.molecular_search_settings.usedAtoms['O'] = (1, 20)
        mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0, 4)
        #mass_spectrum.molecular_search_settings.usedAtoms['S'] = (0, 0)
        #mass_spectrum.molecular_search_settings.usedAtoms['P'] = (0, 1)
        #mass_spectrum.molecular_search_settings.usedAtoms['Na'] = (0, 1)

        mass_spectrum.molecular_search_settings.isProtonated = True
        mass_spectrum.molecular_search_settings.max_oc_filter=1.2
        mass_spectrum.molecular_search_settings.max_hc_filter=3
        mass_spectrum.molecular_search_settings.used_atom_valences = {'C': 4,
                                                                        '13C': 4,
                                                                        'H': 1,
                                                                        'O': 2,
                                                                        'N': 3,
                                                                        'P': 3,
                                                                        'S': 2,
                                                                        'Na': 1}


        SearchMolecularFormulas(mass_spectrum, first_hit=True).run_worker_mass_spectrum()
        mass_spectrum.percentile_assigned(report_error=True)

        assignments=mass_spectrum.to_dataframe()
        assignments['Time']=timestart
        results.append(assignments)

    results=pd.concat(results,ignore_index=True)

    return(results)

def error_plots(output,f):
    fig, ((ax1, ax2)) = plt.subplots(1,2)
    fig.set_size_inches(12, 6)
    sns.scatterplot(x='m/z',y='m/z Error (ppm)',hue='Molecular class',data=output,ax=ax1, edgecolor='none')
    ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
    ax1.set_title('a', fontweight='bold', loc='left')
    sns.kdeplot(x='m/z Error (ppm)',data=output,hue='Time',ax=ax2,legend=False)
    ax2.set_title('b', fontweight='bold', loc='left')
    fig.tight_layout()
    fig_name = f.replace('.raw','_errorplot.jpg')
    fig.savefig(fig_name, dpi=200,format='jpg')



if __name__ == '__main__':

    data_dir = '/CoreMS/usrdata/'
    results = []

    #Set time ranges here
    interval = 2
    time_min = 2
    time_max = 36

    
    #Define time ranges to search
    times = list(range(time_min,time_max,interval))

    #Find all .raw files in the data_dir
    flist = os.listdir(data_dir)
    f_raw = [f for f in flist if '.raw' in f]

    # Generate molecular formula assignments for each file. 
    for f in f_raw:
        print(f)
        output = assign_formula(file = f, times = times)
        output['File'] = f
        output['Molecular class']=output['Molecular Formula'].str.replace('\d+', '').str.replace(' ', '')

        # Save output
        fname = f.replace('.raw','_assigned.csv')
        output.to_csv(data_dir+fname)

        # Plot and save error distribution figure
        error_plots(output,f)