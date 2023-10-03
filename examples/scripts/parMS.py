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

#from corems.mass_spectra.input import rawFileReader
#from corems.molecular_id.factory.classification import HeteroatomsClassification, Labels
#from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
#from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
#from corems.encapsulation.factory.parameters import MSParameters
#from corems.encapsulation.constant import Atoms

from multiprocessing import Pool, Process


import os

num_processor=4

'''

def info(title):
    print(title)
    print('module name:', __name__)
    print('parent process:', os.getppid())
    print('process id:', os.getpid())

def f(name):
    info('function f')
    print('hello', name)

if __name__ == '__main__':
    info('main line')
    p = Process(target=f, args=('bob',))
    p.start()
    p.join()
'''

def get_esi_parser(esifile):
    print(esifile)
    #esi_parser = rawFileReader.ImportMassSpectraThermoMSFileReader(esifile)
    print('module name:', __name__)
    print('parent process:', os.getppid())
    print('process id:', os.getpid())


if __name__ == '__main__':
    esifile = '/Users/christiandewey/CoreMS/tests/tests_data/marine-iodine/220822_CTD27_600m2.raw'
   # get_esi_parser(esifile)
    p = Process(target=get_esi_parser,args=(esifile,))
    p.start()
    p.join()





