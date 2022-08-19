#### Matches LC-ICPMS data and LC-ESIMS data and assignes MF to peak.

#### Christian Dewey
#### 18 August 2022

import os
import pandas as pd
import numpy as np

print(os.path.abspath(os.curdir))
corems_dir = os.path.abspath(os.curdir)
# Change the current working directory
#os.chdir(os.path.dirname(os.getcwd()))
#print(os.path.abspath(os.curdir))

#os.chdir('/Users/boiteaur/Desktop/CoreMS_metallomics/CoreMS')

import warnings
warnings.filterwarnings("ignore")

import sys
sys.path.append("./")

from pathlib import Path

import matplotlib.pyplot as plt
# from PySide2.QtWidgets import QFileDialog, QApplication
# from PySide2.QtCore import Qt

from corems.mass_spectra.input import rawFileReader





from corems.lc_icpms_ftms.factory.Aligned_ICP_ESI import Aligned_ICP_ESI as alignedMS

esif = '/Users/christiandewey/CoreMS/tests/tests_data/ftms/rmb_161221_kansas_h2o_2.raw'
icpf = '/Users/christiandewey/CoreMS/tests/tests_data/icpms/161220_soils_hypercarb_3_kansas_qH2O.csv'

test = alignedMS(icpf,esif)

test.getMSData()
subicp = test.subset_icpdata()

mzs, avms = test.subset_esidata()

results = test.assignFormulas()
