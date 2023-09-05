import os
from tempfile import tempdir
import time
import numpy as np
import warnings
from datetime import date, datetime
import pandas as pd

warnings.filterwarnings('ignore')
from pathlib import Path
import sys
sys.path.append('./')

os.chdir('/Users/christiandewey/CoreMS/')
from corems.mass_spectra.input import rawFileReader
from corems.molecular_id.factory.classification import HeteroatomsClassification, Labels
from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.encapsulation.factory.parameters import MSParameters
from corems.encapsulation.constant import Atoms
from corems.mass_spectrum.calc.Calibration import MzDomainCalibration
os.chdir('/Users/christiandewey/data-temp/tricho/test')


def assign_formula(esifile, times):

	#global search settings
	MSParameters.molecular_search.url_database = 'postgresql+psycopg2://coremsappdb:coremsapppnnl@localhost:5432/coremsapp'
	MSParameters.molecular_search.error_method = None
	MSParameters.molecular_search.min_ppm_error = -1
	MSParameters.molecular_search.max_ppm_error = 1
	MSParameters.molecular_search.score_method = 'prob_score'
	MSParameters.molecular_search.output_score_method = 'prob_score'

	MSParameters.molecular_search.default_ion_charge = 1
	MSParameters.molecular_search.min_ion_charge = 1
	MSParameters.molecular_search.max_ion_charge = 2
	
	MSParameters.mass_spectrum.min_picking_mz=100
	MSParameters.mass_spectrum.max_picking_mz=800
	MSParameters.mass_spectrum.threshold_method = "log"
	MSParameters.mass_spectrum.log_nsigma=400
	MSParameters.ms_peak.peak_min_prominence_percent = 0.1

	# calibration settings
	MSParameters.mass_spectrum.min_calib_ppm_error = -10
	MSParameters.mass_spectrum.max_calib_ppm_error = 10
	MSParameters.mass_spectrum.calib_pol_order = 2
	MSParameters.mass_spectrum.calib_sn_threshold = 5


	parser = rawFileReader.ImportMassSpectraThermoMSFileReader(esifile)

	tic=parser.get_tic(ms_type='MS')[0]
	tic_df=pd.DataFrame({'time': tic.time,'scan': tic.scans})
	results = []
	refmasslist = "/Users/christiandewey/CoreMS/db/nom_pos.ref"
	
	for timestart in times:

		scans=tic_df[tic_df.time.between(timestart,timestart+interval)].scan.tolist()
		mass_spectrum = parser.get_average_mass_spectrum_by_scanlist(scans) 

		#MzDomainCalibration(mass_spectrum, refmasslist,mzsegment=[100,800]).run()

		#first search settings
		mass_spectrum.molecular_search_settings.min_dbe = 0
		mass_spectrum.molecular_search_settings.max_dbe = 30
		mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1, 50)
		mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4, 100)
		mass_spectrum.molecular_search_settings.usedAtoms['O'] = (1, 20)
		mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0, 8)
		mass_spectrum.molecular_search_settings.usedAtoms['S'] = (0, 2)
		mass_spectrum.molecular_search_settings.usedAtoms['P'] = (0, 2)


		mass_spectrum.molecular_search_settings.isProtonated = True
		mass_spectrum.molecular_search_settings.isRadical = False
		mass_spectrum.molecular_search_settings.isAdduct = False
		mass_spectrum.molecular_search_settings.used_atom_valences = {'C': 4,
		                                                                '13C': 4,
		                                                                'H': 1,
		                                                                'O': 2,
		                                                                'N': 3,
                                                                        'S': 2,
                                                                        'P': 3}

		SearchMolecularFormulas(mass_spectrum,first_hit = False).run_worker_mass_spectrum()


		mass_spectrum.percentile_assigned(report_error=True)
		assignments=mass_spectrum.to_dataframe()
		assignments['Time']=timestart
		results.append(assignments)

	results=pd.concat(results,ignore_index=True)

	return(results)



if __name__ == '__main__':

	data_dir = '/Users/christiandewey/data-temp/tricho/test/'
	results = []

	interval = 2
	time_min = 2
	time_max = 30

	times = list(range(time_min,time_max,interval))

	flist = os.listdir(data_dir)
	f_raw = [f for f in flist if '.raw' in f]
	os.chdir(data_dir)
	i=1

	for f in f_raw:
		print(f)
		output = assign_formula(esifile = f, times = times)
		output['file'] = f
		results.append(output)
		i = i + 1 

	fname = 'assignments.csv'
	df = pd.concat(results)
	df.to_csv(data_dir+fname)
