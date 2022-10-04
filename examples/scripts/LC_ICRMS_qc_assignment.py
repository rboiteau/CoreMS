import os
from tempfile import tempdir
from termios import TOSTOP
from time import time
from turtle import color
from unittest.mock import NonCallableMagicMock
import pandas as pd
import numpy as np
import warnings
import math
import re
warnings.filterwarnings("ignore")

import sys
sys.path.append("./")

from pathlib import Path

import tqdm
import importlib

import matplotlib.pyplot as plt
# from PySide2.QtWidgets import QFileDialog, QApplication
# from PySide2.QtCore import Qt
from copy import copy 

from corems.mass_spectra.input import rawFileReader
from corems.molecular_id.factory.classification import HeteroatomsClassification, Labels
from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.encapsulation.factory.parameters import MSParameters
from corems.encapsulation.constant import Atoms


import importlib

import corems.lc_icpms_ftms.calc.lc_icrms_qc_assign as icrms

importlib.reload(icrms)

#Set peak detection threshold method
MSParameters.mass_spectrum.threshold_method = 'signal_noise'
MSParameters.mass_spectrum.s2n_threshold = 2
MSParameters.ms_peak.peak_min_prominence_percent = 0.001

MSParameters.molecular_search.error_method = 'None'
MSParameters.molecular_search.min_ppm_error = -5
MSParameters.molecular_search.max_ppm_error = 5

MSParameters.molecular_search.url_database = None
MSParameters.molecular_search.min_dbe = -1
MSParameters.molecular_search.max_dbe = 20

MSParameters.molecular_search.usedAtoms['C'] = (1,50)
MSParameters.molecular_search.usedAtoms['H'] = (4,100)
MSParameters.molecular_search.usedAtoms['O'] = (1,20)
MSParameters.molecular_search.usedAtoms['N'] = (0,4)
#MSParameters.molecular_search.usedAtoms['S'] = (0,0)

MSParameters.molecular_search.isProtonated = True
MSParameters.molecular_search.isRadical = False
MSParameters.molecular_search.isAdduct = False

MSParameters.molecular_search.score_method = "prob_score"
MSParameters.molecular_search.output_score_method = "prob_score"


# create data class and load parsers 
file_location='/Users/christiandewey/Downloads/Hawaiian_soils/Thermo_RAW/pos/'

hawaii_data = icrms.lc_icr_assign(file_location)


##### get EIC of b12 std; use to view b12 peaks and determine retention time (visual inspection of plot)
hawaii_data.get_b12_eic() 


##### run b12 qc 
hawaii_data.run_b12_qc(b12_peakrange = [8,8.5]) 


#### assign formula 
timerange = [2,10]  # start (min) to stop (min)
interval = 5  # min 
refmasslist = Path.cwd() / "tests/tests_data/ftms/nom_pos.ref"

hawaii_data.assign_formula(interval = interval, timerange = timerange, refmasslist = refmasslist)


# plot req'd resolving power 
hawaii_data.plot_reqd_resolving_power()


# plot results 
hawaii_data.assess_all_results()


# plot unique results 
hawaii_data.determine_unique_features()


# set sort order for clustering analysis
sort_order=[]

for file in hawaii_data._raw_filelist:

    if ('fp_p' in file):  # this portion is unique to naming convention 
        d=file.split('_')[3] 
        sort_order.append(d)


# run cluster analysis
hawaii_data.run_cluster_analysis(sort_order)