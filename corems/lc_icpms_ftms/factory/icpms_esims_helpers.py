from termios import TOSTOP
from unittest.mock import NonCallableMagicMock
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")

import sys
sys.path.append("./")

from pathlib import Path

from corems.mass_spectra.input import rawFileReader
from corems.molecular_id.factory.classification import HeteroatomsClassification, Labels
from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.encapsulation.factory.parameters import MSParameters
from corems.encapsulation.constant import Atoms

def get_esiparser(esifile):
    MSParameters.mass_spectrum.threshold_method = 'signal_noise'
    MSParameters.mass_spectrum.s2n_threshold = 6
    MSParameters.ms_peak.peak_min_prominence_percent = 0.1 #0.1
    return(rawFileReader.ImportMassSpectraThermoMSFileReader(esifile))

def extract(ms):
    mzexp = ms[1]
    esi_parser = ms[0]
    EICdic = {}
    #pbar = tqdm.tqdm(, desc="Getting EICs")
    
    for mz in mzexp:   
        #   print('AverageMS mz:' + str(mz))
        EIC=esi_parser.get_eics(target_mzs=[mz],tic_data={},peak_detection=False,smooth=False)
        EICdic[mz]=EIC[0][mz]

    return EICdic