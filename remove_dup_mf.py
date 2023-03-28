import os
from tempfile import tempdir
import time
import numpy as np
import warnings
from datetime import date, datetime
import pandas as pd

warnings.filterwarnings("ignore")
from pathlib import Path
import sys
sys.path.append("./")

from corems.mass_spectra.input import rawFileReader
from corems.molecular_id.factory.classification import HeteroatomsClassification, Labels
from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.encapsulation.factory.parameters import MSParameters
from corems.encapsulation.constant import Atoms
from corems.mass_spectrum.calc.Calibration import MzDomainCalibration
import corems.lc_icpms_ftms.calc.lc_icrms_qc_assign as icrms
import corems.lc_icpms_ftms.calc.lc_icrms_helpers as lcmsfns

if __name__ == '__main__':
    start = time.time()  #for duration
    startdt = datetime.now()
    

    data_dir = '/home/dewey/Rawfiles/iulia-iodine/'
    fname = 'Organics_mix_02_c.csv'

    os.chdir(data_dir)



    df = pd.read_csv(fname)

    df['abserr'] = np.abs(df['m/z Error (ppm)'])

    df.sort_values(by=['abserr'], inplace=True)

    df.drop_duplicates(subset='Molecular Formula', inplace=True)

    newfname = 'Organics_mix_02_processed.csv'
    df.to_csv(data_dir+newfname)


    
