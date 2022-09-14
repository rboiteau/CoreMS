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

from corems.lc_icpms_ftms.factory.Aligned_ICP_ESI import Aligned_ICP_ESI as alignedMS

esifile = '/Users/christiandewey/CoreMS/tests/tests_data/cobalt-wastewater/rmb_20220627_fp_1uM_stdmix_45.raw'
icpmsfile = '/Users/christiandewey/CoreMS/tests/tests_data/cobalt-wastewater/cwd_220627_3_biozen_1uMstdmix_50uL.csv'

esifile = '/Users/christiandewey/CoreMS/tests/tests_data/cobalt-wastewater/rmb_20220627_fp_CWD_wastewater_43.raw'
icpmsfile = '/Users/christiandewey/CoreMS/tests/tests_data/cobalt-wastewater/cwd_220627_42_biozen_wastewater_50uL.csv'

heteroatom = '59Co'

alignedData = alignedMS(icpmsfile,esifile,heteroatom)

MSParameters.mass_spectrum.threshold_method = 'signal_noise'
MSParameters.mass_spectrum.s2n_threshold = 6
MSParameters.ms_peak.peak_min_prominence_percent = 0.1

esi_parser = rawFileReader.ImportMassSpectraThermoMSFileReader(esifile)

alignedData.timestart = 14  #8.0
alignedData.timestop = 15 #8.5
alignedData.offset = 35.5

alignedData.getMSData()
icpsub = alignedData.subset_icpdata()
icp_subset = icpsub


fig, ax = plt.subplots()
ax.plot(icp_subset['Time '+heteroatom],icp_subset[heteroatom])
plt.show()

mzcorr, mass_spectrum = alignedData.subset_esidata()


elementDict = {'C': (0,70), 'H':(0,90), 'O':(0,15), 'N':(0,15), 'P':(0,1), 'Co':(0,1)} 
threshold = 0.3

#results = alignedData.assignFormulas(elements, threshold)


mass_spectrum.molecular_search_settings.error_method = 'None'
mass_spectrum.molecular_search_settings.min_ppm_error = -10
mass_spectrum.molecular_search_settings.max_ppm_error = 10

mass_spectrum.molecular_search_settings.url_database = None
mass_spectrum.molecular_search_settings.min_dbe = 0
mass_spectrum.molecular_search_settings.max_dbe = 100

elementList = elementDict.keys()

for e in elementList:
    mass_spectrum.molecular_search_settings.usedAtoms[e] = elementDict[e]


mass_spectrum.molecular_search_settings.isProtonated = True
mass_spectrum.molecular_search_settings.isRadical = False
mass_spectrum.molecular_search_settings.isAdduct = False


mass_spectrum.molecular_search_settings.ion_charge = 1


# mass_spectrum.filter_by_max_resolving_power(15, 2)
SearchMolecularFormulas(mass_spectrum, first_hit=False).run_worker_mass_spectrum()

mass_spectrum.molecular_search_settings.score_method = "prob_score"
mass_spectrum.molecular_search_settings.output_score_method = "prob_score"
mass_spectrum.percentile_assigned(report_error=True)

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

mzs_corr = mzcorr
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


match = re.findall(r'[A-Za-z]+|\d+', heteroatom)
heteroAtom = match[1]


results=assignments[assignments['corr']>threshold].filter(['m/z','corr','Peak Height','Confidence Score','Molecular Formula',heteroAtom])

print(results)

bestresults=results[(results[heteroAtom]>=1)]


print(bestresults)
