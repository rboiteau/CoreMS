import os
from tempfile import tempdir
import time
import numpy as np
import warnings
from datetime import date, datetime
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

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
os.chdir('/Users/christiandewey/data-temp/bats')


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

	mz_nform = pd.DataFrame(columns=['timestart','mz_exp','n_candidates'])
	
	for timestart in times:

		scans=tic_df[tic_df.time.between(timestart,timestart+interval)].scan.tolist()
		mass_spectrum = parser.get_average_mass_spectrum_by_scanlist(scans) 

		MzDomainCalibration(mass_spectrum, refmasslist,mzsegment=[100,800]).run()

		#first search settings
		mass_spectrum.molecular_search_settings.min_dbe = 0
		mass_spectrum.molecular_search_settings.max_dbe = 20
		mass_spectrum.molecular_search_settings.error_method = 'None'
		mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1, 40)
		mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4, 80)
		mass_spectrum.molecular_search_settings.usedAtoms['O'] = (0, 20)
		mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0, 8)
		mass_spectrum.molecular_search_settings.usedAtoms['S'] = (0, 2)
		mass_spectrum.molecular_search_settings.usedAtoms['P'] = (0, 2)
		mass_spectrum.molecular_search_settings.usedAtoms['Cu'] = (0, 1)


		mass_spectrum.molecular_search_settings.isProtonated = True
		mass_spectrum.molecular_search_settings.isRadical = False
		mass_spectrum.molecular_search_settings.isAdduct = False
		mass_spectrum.molecular_search_settings.used_atom_valences = {'C': 4,
		                                                                '13C': 4,
		                                                                'H': 1,
		                                                                'O': 2,
		                                                                'N': 3,
                                                                        'S': 2,
                                                                        'P': 3,
																		'Cu':2 }

		SearchMolecularFormulas(mass_spectrum,first_hit = False).run_worker_mass_spectrum()

		for mspeak in mass_spectrum.sort_by_mz():
			mz_nform = mz_nform.append({'timestart': timestart, 'mz_exp':mspeak.mz_exp, 'n_candidates': mspeak.n_candidates},ignore_index = True) 

		mass_spectrum.percentile_assigned(report_error=True)
		assignments=mass_spectrum.to_dataframe()
		assignments['Time']=timestart
		results.append(assignments)

	results=pd.concat(results,ignore_index=True)

	mz_nform.to_csv('n_candidates.csv', index = False)

	_, ax = plt.subplots()

	#ax = sns.barplot(data = mz_nform, x = 'mz_exp', y = 'n_candidates', ax = ax)
	ax.bar(mz_nform['mz_exp'],mz_nform['n_candidates'])
	ax.set_xlim(100,800)
	ax.xaxis.set_major_locator(ticker.MultipleLocator(100))
	
	plt.savefig('barplt.pdf') 
	
	return(results)



if __name__ == '__main__':

	data_dir = '/Users/christiandewey/data-temp/bats/'
	results = []

	interval = 2
	time_min = 10
	time_max = 12

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
