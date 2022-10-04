import os
from termios import TOSTOP
from unittest.mock import NonCallableMagicMock
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import sys
sys.path.append("./")
from pathlib import Path
import matplotlib.pyplot as plt

from corems.mass_spectra.input import rawFileReader
#from corems.molecular_id.factory.classification import HeteroatomsClassification, Labels
#from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
#from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
#from corems.encapsulation.factory.parameters import MSParameters
#from corems.encapsulation.constant import Atoms

from multiprocessing import Pool
import multiprocessing
import time

from tqdm import *
import os



def _run_eics(esifile, mz):

    esi_parser = rawFileReader.ImportMassSpectraThermoMSFileReader(esifile)
    
   # mass_spectrum = esi_parser.get_average_mass_spectrum_by_scanlist(slist)
    
    #rlabel = "%s-%s" %(slist[0],slist[-1])
   # print(rlabel)
    EICdic = {}
    
   # for mz in mz_exp:   
    
    EIC=esi_parser.get_eics(target_mzs=[mz],tic_data={},peak_detection=False,smooth=False)

    EICdic[mz]=EIC[0][mz]

    return EICdic

def _assign_mz(esifile, mz):

    esi_parser = rawFileReader.ImportMassSpectraThermoMSFileReader(esifile)
    
   # mass_spectrum = esi_parser.get_average_mass_spectrum_by_scanlist(slist)
    
    #rlabel = "%s-%s" %(slist[0],slist[-1])
   # print(rlabel)
    EICdic = {}
    
   # for mz in mz_exp:   
    
    EIC=esi_parser.get_eics(target_mzs=[mz],tic_data={},peak_detection=False,smooth=False)

    EICdic[mz]=EIC[0][mz]

    return EICdic

def get_eics(esifile, mz_exp, nprocessors=None):
    
    start_time = time.time()
    
    args_holder = []

    if nprocessors:
        
        nps = nprocessors
    
    else:
        nps = multiprocessing.cpu_count() 

    for mz in mz_exp:

        args_holder.append((esifile,mz))
    
    with Pool(processes=nps) as pool:

        output = pool.starmap(_run_eics, args_holder)

    results = {}
    
    for l in output:
        results = {**results, **l} 
    

    print("--- %s EICs extracted in %.2f seconds using %s cores ---" % (len(output), (time.time() - start_time), nps))

    return results
    
if __name__ == '__main__':
    #esifile = '/Users/christiandewey/CoreMS/tests/tests_data/cobalt-wastewater/rmb_20220627_fp_CWD_wastewater_43.raw'
    
    #slist = [[2255], [2263]] #, 2269] #, [2278, 2288, 2296], [2310, 2319, 2331], [2338, 2343, 2354], [2366, 2374, 2383], [2398, 2411, 2418], [2422, 2433, 2441], [2448, 2454, 2466]]

    #args_holder = []
    #for ll in slist:
    #    args_holder.append((esifile,ll))
    num_processors = 8

    #with Pool(processes=num_processors) as pool:
    #    results = pool.starmap(get_eics, args_holder)
    

