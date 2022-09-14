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


os.chdir('tests/tests_data/marine-iodine')
os.getcwd()

fname = '220822_CTD27_1000m2.raw'

MSParameters.mass_spectrum.threshold_method = 'signal_noise'
MSParameters.mass_spectrum.s2n_threshold = 6
MSParameters.ms_peak.peak_min_prominence_percent = 0.1
esi_parser =  rawFileReader.ImportMassSpectraThermoMSFileReader(fname)


# peak window for b12 

t1 = 34.75# 35.24 # min
t2 = 37.25 #36.83 # min

tic = esi_parser.get_tic(ms_type='MS')[0]


tic_df = pd.DataFrame({'time': tic.time, 'scan':tic.scans})
scans = tic_df[tic_df.time.between(t1,t2)].scan.tolist()

averagems = esi_parser.get_average_mass_spectrum_by_scanlist(scans)

eicdic = {}
pbar = tqdm.tqdm(averagems.mz_exp)
for mz in pbar:
    eic = esi_parser.get_eics(target_mzs=[mz],tic_data={},peak_detection=False,smooth=False)
    eicdic[mz] = eic[0][mz]
    

for key in eicdic.keys():
    print(key)

fig, ax = plt.subplots()

x = eicdic[678.282385619928].time
y = eicdic[678.282385619928].eic
ax.plot(x,y)

ax.set_xlim((35,38))

plt.show()

times=tic_df[tic_df.time.between(t1,t2)].time.tolist()


### 35.85 --> icpms b12 peak
### 35.20 --> icrms b12 peak

## icrms is 39 s ahead of icpms 
