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
    ax.text(0,m, 'ratio: ' + ratio)
    ax.text(0,m*.9, 'charge: ' + str(candidate.charge))
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




def isotope_hunter(ms, tolerance, dominant, secondary):
    m1 = Atoms.atomic_masses[secondary]# ['Fe'] #
    m2 = Atoms.atomic_masses[dominant]  # ['54Fe'] #
    #m3 = Atoms.atomic_masses['57Fe']
    #m4 = Atoms.atomic_masses['58Fe']

    candidates = [] 
    holder = []
    
   # for charge in charge_range:
    diff1 = abs(m2 - m1 )
    diff2 = abs(m2 - m1) / 2 
    diff3 = abs(m2 - m1) / 3

   # print('Mass difference between ' + dominant + ' and ' + secondary + ' (charge = +1): ' + str(diff1))
   # print('Mass difference between ' + dominant + ' and ' + secondary + ' (charge = +2): ' + str(diff2))
    #print('Mass difference between ' + heavy + ' and ' + light + ' (charge = +3): ' + str(diff3))
    carbon_diff = Atoms.atomic_masses['13C'] - Atoms.atomic_masses['C']

   #     print('charge: ' + str(charge))

    for mz1, i in zip(ms.mz_exp, range(len(ms.mz_exp))):

    
        #if (mz1 < 384) and (mz1 > 380):
        #    print('mz1 %.4f' % mz1)
        for mz2, a in zip(ms.mz_exp[i+1:], ms.abundance[i+1:]):
            
            if m2 > m1:

                c_dom = mz2

                j = i + 2 

            elif m1 > m2:

                c_dom = mz1

                j = i 

            if abs(mz2 - mz1) > (diff1 * 1.1):
                break 

            elif (abs(mz2 - mz1) < (tolerance + diff1)) and (abs(mz2 - mz1)  > ( (-1 * tolerance) + diff1)):      
          #      hg = abs(mz2-mz1)
          #      print('\tabs(mz2-mz1): %.4f' % hg)
          #      hh = tolerance + diff1
          #      print('\ttolerance + diff1: %.4f' % hh)
          #      hj = diff1 - tolerance
          #      print('\t-tolerance + diff1: %.4f' % hj)

                try: 
                
                    while ms.mz_exp[j] < (c_dom + (carbon_diff * 1.1) ):
                   #     print('\t\tmz[j-1]  ' + str(abs(ms.mz_exp[j-1])))
                    #    print('\t\tmz[j]  ' + str(abs(ms.mz_exp[j])))
                     #   print('\t\tmz[j+1]  ' + str(abs(ms.mz_exp[j+1])))
                      #  print('\t\tmz[j+2]  ' + str(abs(ms.mz_exp[j+2])))
                       # gh = tolerance +  carbon_diff  
                       # print('\t\ttolerance +  carbon_diff: %.4f' % gh)
                       # gj = abs(ms.mz_exp[j] - c_dom)
                       # print('\t\tabs(ms.mz_exp[j] - c_dom): %.4f' % gj)
                        
                        if ( abs(ms.mz_exp[j] - c_dom) < ( tolerance +  carbon_diff   ) ) and ( abs(ms.mz_exp[j] - c_dom) > ( ( -1 * tolerance)  +  carbon_diff ) ):

                            print('charge: 1\t' + 'm2: ' + str(mz2) + '\tm1: ' + str(mz1) + '\tdiff: ' + str(mz2-mz1) + '\ta1: ' + str(ms.abundance[i]) + '\ta2: ' + str(a))

                            candidates.append(mz2) 

                            charge = 1

                            clist = [mz1, mz2]

                            alist = [ms.abundance[i], a]

                            eiclist = [esi_parser.get_eics(target_mzs=[mz1],tic_data={},peak_detection=False,smooth=False), esi_parser.get_eics(target_mzs=[mz2],tic_data={},peak_detection=False,smooth=False)]

                            holder.append(Candidate(clist,alist,eiclist,charge))
                        
                        j = j + 1

                except:
                    
                    pass

            elif (abs(mz2 - mz1) < (tolerance + diff2)) and (abs(mz2 - mz1)  > ( (-1 * tolerance) + diff2)):

                try:
                    
                    while ms.mz_exp[j] < (c_dom + ( carbon_diff / 2 ) * 1.1):
                        
                        if ( abs(ms.mz_exp[j] - c_dom) < ( tolerance + ( carbon_diff / 2 ) ) ) and ( abs(ms.mz_exp[j] - c_dom) > ( ( -1 * tolerance)  + ( carbon_diff / 2 ) ) ) :

                            print('charge: 2\t' + 'm2: ' + str(mz2) + '\tm1: ' + str(mz1) + '\tdiff: ' + str(mz2-mz1) + '\ta1: ' + str(ms.abundance[i]) + '\ta2: ' + str(a))

                            candidates.append(mz2) 

                            charge = 2

                            clist = [mz1, mz2]

                            alist = [ms.abundance[i], a]

                            eiclist = [esi_parser.get_eics(target_mzs=[mz1],tic_data={},peak_detection=False,smooth=False), esi_parser.get_eics(target_mzs=[mz2],tic_data={},peak_detection=False,smooth=False)]

                            holder.append(Candidate(clist,alist,eiclist,charge))
                        
                        j = j + 1
                
                except:

                    pass

            elif (abs(mz2 - mz1) < (tolerance + diff3)) and (abs(mz2 - mz1)  > ( (-1 *tolerance) + diff3)):

                try: 

                    while ms.mz_exp[j] < (c_dom + ( carbon_diff / 3 )+ 1.1*tolerance):
                        
                        if ( abs(ms.mz_exp[j] - c_dom) <( tolerance + ( carbon_diff / 3 ) ) ) and ( abs(ms.mz_exp[j] - c_dom) > ( ( -1 * tolerance)  + ( carbon_diff / 3 ) ) ):
                            
                            print('charge: 3\t' + 'm2: ' + str(mz2) + '\tm1: ' + str(mz1) + '\tdiff: ' + str(mz2-mz1) + '\ta1: ' + str(ms.abundance[i]) + '\ta2: ' + str(a))

                            candidates.append(mz2) 

                            charge = 3

                            clist = [mz1, mz2]

                            alist = [ms.abundance[i], a]

                            eiclist = [esi_parser.get_eics(target_mzs=[mz1],tic_data={},peak_detection=False,smooth=False), esi_parser.get_eics(target_mzs=[mz2],tic_data={},peak_detection=False,smooth=False)]

                            holder.append(Candidate(clist,alist,eiclist,charge))

                        j = j + 1

                except:
                    
                    pass

        dict = {}
        
        for c, a in zip(candidates, holder):
            
            dict[c] = a

    return dict 



esifile = '/Users/christiandewey/Downloads/rmb_20220627_fp_1uM_stdmix_45.raw'
#esifile = '/Users/christiandewey/CoreMS/tests/tests_data/marine-iodine/220822_CTD27_600m.raw'
#esifile ='/Users/christiandewey/CoreMS/tests/tests_data/cobalt-wastewater/rmb_20220627_fp_CWD_wastewater_43.raw'
#esifile = '/Users/christiandewey/CoreMS/tests/tests_data/ftms/rmb_161221_kansas_h2o_2.raw'


MSParameters.mass_spectrum.threshold_method = 'signal_noise'
MSParameters.mass_spectrum.s2n_threshold = 3# 6
MSParameters.ms_peak.peak_min_prominence_percent = 0.1 #0.1

esi_parser = rawFileReader.ImportMassSpectraThermoMSFileReader(esifile)

tic = esi_parser.get_tic(ms_type = 'MS')[0]

#plot_eic(tic.time, tic.tic)

#plt.show()#


scan_index = np.column_stack((range(len(tic.scans)), tic.scans, tic.time))

tstart = 0 #min
tstop = max(tic.time)


chroma_window = 20 # sec

scan_freq  = tic.time[10] - tic.time[9] # min

scans_per_chroma = int ( chroma_window / (scan_freq * 60) ) + 1

sub = scan_index[np.where((scan_index[:,2]>tstart) & (scan_index[:,2]<=tstop))]

for scan, time in zip(sub[:,1], sub[:,2]):
    print(scan, time)

temp = list(sub[:,1])
scanlist = [int(i) for i in temp]

scanlist
srange = range(0,len(scanlist),scans_per_chroma)
for r in srange:
    print(r)
candidates = []
all_mz = []

for r in srange: 
    
    scanrange = scanlist[r:r+scans_per_chroma]

    print(scanrange)
   
    scan_ms = esi_parser.get_average_mass_spectrum_by_scanlist(scanrange)
   # all_mz = all_mz + list(scan_ms.mz_exp)

    #plot_ms(scan_ms)

    #plt.show()

    d1 = isotope_hunter(scan_ms, tolerance=0.0005, dominant='Fe', secondary='54Fe')

    if d1.keys():
        candidates.append(d1)


all_mz.sort()




for c in candidates:

    for k, v in zip(c.keys(),c.values()):

        print(v.charge, '%.4f' % v.mz1, '%.4f' % v.mz2)

for c in candidates:

    for k, v in zip(c.keys(),c.values()):

        if v.ratio > 1:

            plot_candidate(v)

            plt.show()





