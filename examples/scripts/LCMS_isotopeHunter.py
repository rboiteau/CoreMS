import os
from termios import TOSTOP
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
from corems.encapsulation.constant import Atoms



class Candidate:

    def __init__(self, clist, alist, eiclist, charge): #,esi_file_location, heteroatom):

        self.mz1 = clist[0]
        self.mz2 = clist[1]

        self.ab1 = alist[0]
        self.ab2 = alist[1]

        self.eic1 = eiclist[0]
        self.eic2 = eiclist[1]

        self.charge = charge

        self._getRatio()

    def _getRatio(self):

        self.ratio = (self.ab2 / self.ab1) / (Atoms.isotopic_abundance['Fe']/Atoms.isotopic_abundance['54Fe'])




def plot_candidate(candidate):

    ax = plt.gca()

    ax.plot(candidate.eic1[0][candidate.mz1].time, candidate.eic1[0][candidate.mz1].eic,label = '%.4f' % candidate.mz1)
    ax.plot(candidate.eic2[0][candidate.mz2].time, candidate.eic2[0][candidate.mz2].eic,label = '%.4f' % candidate.mz2)
    ratio = '%.4f' % candidate.ratio
    m = max(candidate.eic1[0][candidate.mz1].eic + candidate.eic2[0][candidate.mz2].eic)
    ax.text(0,m, ratio)
    ax.text(0,m*.9, candidate.charge)
    ax.legend()

    return ax


def plot_ms(ms):
    mz_exp = ms.mz_exp
    abundance = ms.abundance
    ax = plt.gca()

    for plot_obj in ax.stem(mz_exp, abundance, linefmt = '-', markerfmt=" ", use_line_collection = True, label = 'm/z'):
        plt.setp(plot_obj, 'linewidth', 2)

    return ax



def plot_eic(time, eic):

    ax = plt.gca()

    ax.plot(time, eic)

    return ax 




def isotope_hunter(ms, tolerance):
    m1 = Atoms.atomic_masses['54Fe']# ['Fe'] #
    m2 = Atoms.atomic_masses['Fe']  # ['54Fe'] #
    m3 = Atoms.atomic_masses['57Fe']
    m4 = Atoms.atomic_masses['58Fe']

    i1 = Atoms.isotopic_abundance['Fe']
    i2 = Atoms.isotopic_abundance['54Fe']
    i3 = Atoms.isotopic_abundance['57Fe']
    i4 = Atoms.isotopic_abundance['58Fe']

    em = Atoms.electron_mass
    
    charge_range = [1,2,3]
    candidates = [] 
    holder = []
    
    for charge in charge_range:
        diff = m2 - m1 - ( em * charge )
        

        print('charge: ' + str(charge))

        for mz1, i in zip(ms.mz_exp, range(len(ms.mz_exp))):

            for mz2, a in zip(ms.mz_exp[i+1:], ms.abundance[i+1:]):

                if (abs(mz2 - mz1 - diff) <= tolerance):

                    #candidates.append(mz1)

                    candidates.append(mz2)

                    #abundances.append(ms.abundance[i+1])

        #            abundances.append(a)

                    clist = [mz1, mz2]

                    alist = [ms.abundance[i+1], a]

                    eiclist = [esi_parser.get_eics(target_mzs=[mz1],tic_data={},peak_detection=False,smooth=False), esi_parser.get_eics(target_mzs=[mz2],tic_data={},peak_detection=False,smooth=False)]

                    holder.append(Candidate(clist,alist,eiclist,charge))
            
                elif abs(mz1 - mz2) > (diff + tolerance):

                    break

        dict = {}
        
        for c, a in zip(candidates, holder):
            
            dict[c] = a

    return dict 



esifile = '/Users/christiandewey/Downloads/rmb_20220627_fp_1uM_stdmix_45.raw'
#esifile = '/Users/christiandewey/CoreMS/tests/tests_data/marine-iodine/220822_CTD27_600m.raw'

MSParameters.mass_spectrum.threshold_method = 'signal_noise'
MSParameters.mass_spectrum.s2n_threshold = 6
MSParameters.ms_peak.peak_min_prominence_percent = 0.1

esi_parser = rawFileReader.ImportMassSpectraThermoMSFileReader(esifile)

tic = esi_parser.get_tic(ms_type = 'MS')[0]



scan_index = np.column_stack((range(len(tic.scans)), tic.scans, tic.time))

tstart = 0 #min
tstop = max(tic.time)

sub = scan_index[np.where((scan_index[:,2]>tstart) & (scan_index[:,2]<=tstop))]

for scan, time in zip(sub[:,1], sub[:,2]):
    print(scan, time)

temp = list(sub[:,1])
scanlist = [int(i) for i in temp]

scan_ms = esi_parser.get_average_mass_spectrum_by_scanlist(scanlist)



d1 = isotope_hunter(scan_ms, tolerance=0.0005)

c2 = d1.keys()


for k, v in zip(d1.keys(),d1.values()):

    if v.ratio > 1:

        plot_candidate(v)

        plt.show()



