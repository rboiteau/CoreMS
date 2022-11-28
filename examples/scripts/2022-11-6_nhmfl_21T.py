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

import seaborn as sns
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
from corems.mass_spectrum.calc.Calibration import MzDomainCalibration



import importlib


import corems.lc_icpms_ftms.calc.lc_icrms_qc_assign as icrms

importlib.reload(icrms)
for t in df['Time']:
    print(t)




def getParser(file):
    parser = rawFileReader.ImportMassSpectraThermoMSFileReader(file)
    return parser

def assign_formula(parser, interval, timerange, refmasslist=None):
    #Function to build formula assignment lists
    #Retrieve TIC for MS1 scans over the time range between 'timestart' and 'timestop' 

    tic=parser.get_tic(ms_type='MS')[0]
    tic_df=pd.DataFrame({'time': tic.time,'scan': tic.scans})

    times=list(range(timerange[0],timerange[1],interval))

    results=[]
    
    for timestart in times:

        scans=tic_df[tic_df.time.between(timestart,timestart+interval)].scan.tolist()

        mass_spectrum = parser.get_average_mass_spectrum_by_scanlist(scans)    
        mass_spectrum.molecular_search_settings.ion_charge = 1

        ppp = mass_spectrum.polarity
        print('polarity: %s' %ppp)
        #mass_spectrum.mass_spectrum.settings.calib_sn_threshold
        #mass_spectrum.mass_spectrum.settings.calib_pol_order
        #mass_spectrum.recalibrate_mass_spectrum(mass_spectrum, imzmeas, mzrefs, order=2)
        #MzDomainCalibration(mass_spectrum, ref_file_location).run()

        if refmasslist:
            mass_spectrum.settings.min_calib_ppm_error = 10
            mass_spectrum.settings.max_calib_ppm_error = -10
            calfn = MzDomainCalibration(mass_spectrum, refmasslist)
            ref_mass_list_fmt = calfn.load_ref_mass_list(refmasslist)

            imzmeas, mzrefs = calfn.find_calibration_points(mass_spectrum, ref_mass_list_fmt,
                                                        calib_ppm_error_threshold=(0, 2.0),
                                                        calib_snr_threshold=3)

            calfn.recalibrate_mass_spectrum(mass_spectrum, imzmeas, mzrefs, order=2)


        SearchMolecularFormulas(mass_spectrum, first_hit=False).run_worker_mass_spectrum()

        mass_spectrum.percentile_assigned(report_error=True)

        assignments=mass_spectrum.to_dataframe()

        assignments['Time']=timestart

        results.append(assignments)
    
    results=pd.concat(results,ignore_index=True)

    return(results)    

def filterMzRange(results, mz_range):

    mz_i = mz_range[0]
    mz_f = mz_range[1]

    sub = results[(results['m/z'] >= mz_i) & (results1['m/z'] <= mz_f)]

    return sub

def compareFeature(results1, results2):
    # results 1: full scan
    # results 2: narrow scan
    numrows = len(results2['m/z'])

    temp1 = results1
    temp2 = results2

    temp1['m/z_round'] = np.around(temp1['m/z'],4)
    temp2['m/z_round'] = np.around(temp2['m/z'],4)

    for mz in temp2['m/z_round']:
        if mz in temp1['m/z_round']:
            print(mz)

def plot_ms(masspectrum,srange,target_mz, mzrange,assignment= None, ax_ms=None):   
    mass_spectrum = masspectrum
    if ax_ms == None:
        f, ax = plt.subplots()
    
    else:
        ax = ax_ms

    timerange = srange[1]

    ms_df_all=mass_spectrum.to_dataframe()
    print(ms_df_all)

    ms_df = ms_df_all[(abs(ms_df_all['m/z']-target_mz)<mzrange)]
    print(ms_df)
    mztarget_line = ms_df[ms_df['m/z'] == target_mz]

    print(mztarget_line)    
    _, stemlines1, _ =ax.stem('m/z','Peak Height',data=ms_df,  markerfmt=' ', basefmt=' ', linefmt='blue')
    _, stemlines2, _ =ax.stem('m/z','Peak Height',data=mztarget_line,  markerfmt=' ', basefmt=' ', linefmt='r')
    
    ax.set_ylim(0, max(ms_df['Peak Height']) * 1.1)
    ax.set_xlim(left = target_mz - mzrange, right = target_mz + mzrange) 

    for mzr,peakr in zip(ms_df['m/z'], ms_df['Peak Height']):

        if (mzr- target_mz)  == 0:
            mz_text = ' m/z\n%.4f' % (mzr)
        else:
            mz_text = r'$\Delta$' + ' m/z\n%.4f' % (mzr- target_mz)
        ax.text(mzr, peakr + 0.02 *max(ms_df['Peak Height']), mz_text, ha = 'center', fontsize = 'xx-small', weight = 'bold')

   # theor_mz=pattern.mdiff+result['mass']
   # theor_int=pattern.ratio*result['abundance']
   # ax.stem(theor_mz,theor_int, basefmt=' ',linefmt='gray')

   # for isotope in pattern.isotope[pattern.requirement=='Y']:
   #     ax.stem('mz','intense',data=result[isotope],  markerfmt=' ', basefmt=' ',linefmt='red')
    if ax == None:

        ax.legend(('other', 'target'),bbox_to_anchor=(1.05, 1.0), loc='upper left',frameon=False)

        if(assignment):

            mf = assignment[0]
            score = assignment[1]
            er = assignment[2]

            ax.text(1.05,0.7,mf,transform=ax.transAxes)
            ax.text(1.05,0.6,'Error (ppm) = %.3f ' %er ,transform=ax.transAxes)
            ax.text(1.05,0.5,'Score = %.3f' %score ,transform=ax.transAxes)

    ax.set(xlabel='m/z',ylabel='Intensity')
    ax.set_title('%.2f' %timerange[0] + ' to %.2f' %timerange[1] +' min', fontsize = 'medium')

    ax.axhline(y=0.0, color='black')
    plt.setp(stemlines1, 'linewidth', 0.75)
    plt.setp(stemlines2, 'linewidth', 0.75)
    plt.tight_layout()
    if ax_ms == None:
        return ax

def plot_mz_error(results_df):
    filtered_results=results_df[(results_df['m/z']<800) & (results_df['S/N']>3)]

    filtered_results['N']=filtered_results['N'].fillna(0)
    filtered_results['O']=filtered_results['O'].fillna(0)
    filtered_results['S']=filtered_results['S'].fillna(0)
    filtered_results['P']=filtered_results['P'].fillna(0)
    #filtered_results['Na']=filtered_results['Na'].fillna(0)

    filtered_results['mol_class']='Unassigned'
    filtered_results['mol_class'][filtered_results['C']>0]='CHO'
    filtered_results['mol_class'][(filtered_results['C']>0) & (filtered_results['N']>0.5)]='CHON'
    filtered_results['mol_class'][(filtered_results['C']>0) & (filtered_results['S']>0.5)]='CHOS'
    filtered_results['mol_class'][(filtered_results['C']>0) & (filtered_results['P']>0.5)]='CHOP'
    #filtered_results['mol_class'][(filtered_results['C']>0) & (filtered_results['Na']>0.5)]='CHONa'
    filtered_results['mol_class'][(filtered_results['C']>0) & (filtered_results['S']>0.5) & (filtered_results['N']>0.5)]='CHONS'
    filtered_results['mol_class'][(filtered_results['C']>0) & (filtered_results['P']>0.5) & (filtered_results['N']>0.5)]='CHONP'
    #filtered_results['mol_class'][(filtered_results['C']>0) & (filtered_results['Na']>0.5) & (filtered_results['N']>0.5)]='CHONNa'
    #filtered_results['mol_class'][(filtered_results['C']>0) & (filtered_results['P']>0.5) & (filtered_results['Na']>0.5) & (filtered_results['N']>0.5)]='CHONPNa'


    results= filtered_results[filtered_results['Is Isotopologue']==0]
    results['O/C']=results['O']/results['C']
    results['H/C']=results['H']/results['C']
    results['N/C']=results['N']/results['C']

def assess_all_results(results_df):

    filtered_results=results_df[(results_df['m/z']<800) & (results_df['S/N']>3)]

    filtered_results['N']=filtered_results['N'].fillna(0)
    filtered_results['O']=filtered_results['O'].fillna(0)
    filtered_results['S']=filtered_results['S'].fillna(0)
    filtered_results['P']=filtered_results['P'].fillna(0)
    #filtered_results['Na']=filtered_results['Na'].fillna(0)

    filtered_results['mol_class']='Unassigned'
    filtered_results['mol_class'][filtered_results['C']>0]='CHO'
    filtered_results['mol_class'][(filtered_results['C']>0) & (filtered_results['N']>0.5)]='CHON'
    filtered_results['mol_class'][(filtered_results['C']>0) & (filtered_results['S']>0.5)]='CHOS'
    filtered_results['mol_class'][(filtered_results['C']>0) & (filtered_results['P']>0.5)]='CHOP'
    #filtered_results['mol_class'][(filtered_results['C']>0) & (filtered_results['Na']>0.5)]='CHONa'
    filtered_results['mol_class'][(filtered_results['C']>0) & (filtered_results['S']>0.5) & (filtered_results['N']>0.5)]='CHONS'
    filtered_results['mol_class'][(filtered_results['C']>0) & (filtered_results['P']>0.5) & (filtered_results['N']>0.5)]='CHONP'
    #filtered_results['mol_class'][(filtered_results['C']>0) & (filtered_results['Na']>0.5) & (filtered_results['N']>0.5)]='CHONNa'
    #filtered_results['mol_class'][(filtered_results['C']>0) & (filtered_results['P']>0.5) & (filtered_results['Na']>0.5) & (filtered_results['N']>0.5)]='CHONPNa'


    results= filtered_results[filtered_results['Is Isotopologue']==0]
    results['O/C']=results['O']/results['C']
    results['H/C']=results['H']/results['C']
    results['N/C']=results['N']/results['C']


    colors = {'CHO':'red', 'CHON':'green', 'CHOS':'blue', 'CHONS':'yellow', 'CHONP':'black', 'CHONNa':'cyan','CHONPNa':'pink','CHONa':'aqua','CHOP':'gray'}

    grouped=results.groupby('mol_class')
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,figsize=(10,6))
    ax1 = sns.violinplot(x="Time", y="O/C", hue='mol_class', data=results, ax=ax1)
    ax1.set(xlabel='Time (min)')
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    ax2 = sns.scatterplot(x='m/z',y='m/z Error (ppm)',hue='mol_class',data=results,ax=ax2)
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    #for key, group in grouped:
    #    group.plot(ax=ax,kind='scatter',x='m/z',y='m/z Error (ppm)',color=colors[key],label=key)

    for key, group in grouped:
        group.plot(ax=ax3,kind='scatter',x='m/z',y='Resolving Power',color=colors[key],label=key)
    ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    assign_summary=[]
    for time in filtered_results['Time'].unique():
        current={}
        current['Time']=time
        for mol_class in sorted(filtered_results['mol_class'].unique()):
            current[mol_class]=len(filtered_results[(filtered_results['mol_class']==mol_class) & (filtered_results['Time']==time)])
        assign_summary.append(current)
    df=pd.DataFrame(assign_summary)

    ax4 = df.plot.bar(x='Time',y=df.columns[1:],stacked=True,ylabel='Peaks',ax=ax4)
    ax4.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    #fig.savefig(self.savedir + 'all_results.pdf')
    fig.tight_layout()
    plt.show()





data1='/Users/christiandewey/Downloads/Slough/20221102_LBA_Boiteau_Zorbax3p5_slough_fullmz_01.raw'
data2='/Users/christiandewey/Downloads/Slough/20221102_LBA_Boiteau_Zorbax3p5_slough_400_500_02.raw'
qh2o='/Users/christiandewey/Downloads/slough/20221102_LBA_Boiteau_Zorbax3p5_qh2o_400_500_03.raw'


refmasslist = '/Users/christiandewey/CoreMS/tests/tests_data/ftms/nom_pos.ref'

savedir = '~/Downloads/'
results1 = pd.read_csv(savedir + 'slough_fullmz.csv')
results2 = pd.read_csv(savedir + 'slough_400-500.csv')
qh2o_results = pd.read_csv(savedir + 'qh2o.csv')


# 1: assign formula to full scan 
MSParameters.mass_spectrum.threshold_method = 'signal_noise'
MSParameters.mass_spectrum.s2n_threshold = 2
MSParameters.ms_peak.peak_min_prominence_percent = 0.001

MSParameters.molecular_search.error_method = 'None'
MSParameters.molecular_search.min_ppm_error = -0.4
MSParameters.molecular_search.max_ppm_error = 0.4

MSParameters.molecular_search.isProtonated = True
MSParameters.molecular_search.isRadical = False
MSParameters.molecular_search.isAdduct = False

MSParameters.molecular_search.score_method = "prob_score"
MSParameters.molecular_search.output_score_method = "prob_score"


MSParameters.molecular_search.url_database = None
MSParameters.molecular_search.min_dbe = -1
MSParameters.molecular_search.max_dbe = 20

MSParameters.molecular_search.usedAtoms['C'] = (1,50)
MSParameters.molecular_search.usedAtoms['H'] = (4,100)
MSParameters.molecular_search.usedAtoms['O'] = (1,20)
MSParameters.molecular_search.usedAtoms['N'] = (0,8)
MSParameters.molecular_search.usedAtoms['S'] = (0,1)
MSParameters.molecular_search.usedAtoms['P'] = (0,1)
#MSParameters.molecular_search.usedAtoms['Na'] = (0,0)





parser1 = getParser(data1)
results1 = assign_formula(parser1,interval=2,timerange=[4,24], refmasslist=refmasslist)

parser2 = getParser(data2)
results2 = assign_formula(parser2,interval=2,timerange=[4,24], refmasslist=refmasslist)

qh2o_parser = getParser(qh2o)
qh2o_results = assign_formula(qh2o_parser,interval=2,timerange=[4,24], refmasslist=refmasslist )


## filter rounded calibrated m/z
results1['m/z 4d'] = np.around(results1['Calibrated m/z'], 4)
results1.to_csv(savedir + 'slough_fullmz.csv')
np.shape(results1)


results2['m/z 4d'] = np.around(results2['Calibrated m/z'], 4)
results2.to_csv(savedir + 'slough_400-500.csv')
np.shape(results2)

qh2o_results['m/z 4d'] = np.around(qh2o_results['Calibrated m/z'], 4)
qh2o_results.to_csv(savedir + 'qh2o.csv')
np.shape(qh2o_results)

sub = filterMzRange(results1,[400,500])
sub['m/z 4d'] = np.around(sub['Calibrated m/z'], 4)
sub.to_csv(savedir + 'slough_fullmz_400-500.csv')
np.shape(sub)

blank_subtracted_n = results2[~results2['m/z 4d'].isin(qh2o_results['m/z 4d'])]
blank_subtracted_n.to_csv(savedir + 'slough_400-500mz_blank_subtract.csv')
np.shape(blank_subtracted_n)

blank_subtracted_full = sub[~sub['m/z 4d'].isin(qh2o_results['m/z 4d'])]
blank_subtracted_full.to_csv(savedir + 'slough_fullmz_blank_subtract.csv')
np.shape(blank_subtracted_full)

overlap = blank_subtracted_n[blank_subtracted_n['m/z 4d'].isin(blank_subtracted_full['m/z 4d'])]
overlap.to_csv(savedir + 'slough_fullmz_400-500_overlap.csv')
np.shape(overlap)

narrow_unique = blank_subtracted_n[~blank_subtracted_n['m/z 4d'].isin(blank_subtracted_full['m/z 4d'])]
narrow_unique.to_csv(savedir + 'slough_400-500_unique.csv')
np.shape(narrow_unique)


## filter molform 
sub = filterMzRange(results1,[400,500])
sub.to_csv(savedir + 'slough_fullmz_400-500.csv')
np.shape(sub)

blank_subtracted_n = results2[~results2['Molecular Formula'].isin(qh2o_results['Molecular Formula'])]
blank_subtracted_n.to_csv(savedir + 'slough_400-500mz_blank_subtract.csv')
np.shape(blank_subtracted_n)

blank_subtracted_full = sub[~sub['Molecular Formula'].isin(qh2o_results['Molecular Formula'])]
blank_subtracted_full.to_csv(savedir + 'slough_fullmz_blank_subtract.csv')
np.shape(blank_subtracted_full)

overlap = blank_subtracted_n[blank_subtracted_n['Molecular Formula'].isin(blank_subtracted_full['Molecular Formula'])]
overlap.to_csv(savedir + 'slough_fullmz_400-500_overlap.csv')
np.shape(overlap)

narrow_unique = blank_subtracted_n[~blank_subtracted_n['Molecular Formula'].isin(blank_subtracted_full['Molecular Formula'])]
narrow_unique.to_csv(savedir + 'slough_400-500_unique.csv')
np.shape(narrow_unique)





## plotting
xlow = 150
xhigh = 2000
ymax = 100000

df = narrow_unique
fig, ax = plt.subplots()
ax.stem(df['m/z'],df['Peak Height'], linefmt='-',  markerfmt=" ", use_line_collection =True)
ax.set_xlim(xlow, xhigh)
ax.set_ylim(0,ymax)
ax.set_title('slough, features unique to 400-500 m/z window')
#ax.legend(frameon=False)
fig.tight_layout()

plt.show()

plt.savefig('slough_400-500_unique.pdf')

df = blank_subtracted_full
fig, ax = plt.subplots()
ax.stem(df['m/z'],df['Peak Height'], linefmt='-',  markerfmt=" ", use_line_collection =True)
ax.set_xlim(xlow, xhigh)
ax.set_ylim(0,ymax)
ax.set_title('slough, full scan (150-2000 m/z)')
#ax.legend(frameon=False)
fig.tight_layout()
plt.savefig('slough_fullmz.pdf')

df = overlap
fig, ax = plt.subplots()
ax.stem(df['m/z'],df['Peak Height'], linefmt='-',  markerfmt=" ", use_line_collection =True)
ax.set_xlim(xlow, xhigh)
ax.set_ylim(0,ymax)
ax.set_title('slough, full/400-500 m/z overlap')
#ax.legend(frameon=False)
fig.tight_layout()
plt.savefig('slough_fullmz_400-500_overlap.pdf')

plt.close('all')

#I space
xlow = 448
xhigh = 450
ymax = 20000
df = blank_subtracted_full
df2 = narrow_unique
fig, ax = plt.subplots()

markerline1, stemlines1, _ = ax.stem(df['m/z'],df['Peak Height'], linefmt='-',  markerfmt=" ",use_line_collection =True, label='150-2000 m/z')
plt.setp(stemlines1, 'color', 'C0', 'linewidth', 2)
markerline2, stemlines2, _ = ax.stem(df2['m/z'],df2['Peak Height'], linefmt='-',  markerfmt=" ",use_line_collection =True, label='unique to 400-500 m/z')
plt.setp(stemlines2, 'color', 'C1', 'linewidth', 2)
ax.set_xlim(xlow, xhigh)
ax.set_ylim(0,ymax)
plt.show()


# S/N space
xlow = 400 #448
xhigh = 500 # 450
ymax = 1500
df = blank_subtracted_full
df2 = blank_subtracted_n

fig, ax = plt.subplots()
markerline1, stemlines1, _ = ax.stem(df['m/z'],df['S/N'], linefmt='-',  markerfmt=" ",use_line_collection =True, label='150-2000 m/z')
plt.setp(stemlines1, 'color', 'C0', 'linewidth', 2)
#markerline2, stemlines2, _ = ax.stem(df2['m/z'],df2['S/N'], linefmt='-',  markerfmt=" ",use_line_collection =True, label='unique to 400-500 m/z')
#plt.setp(stemlines2, 'color', 'C1', 'linewidth', 2)
ax.set_xlim(xlow, xhigh)
ax.set_ylim(0,ymax)
ax.set_xlabel('m/z')
ax.set_ylabel('S/N')
ax.set_title('150-2000 m/z scan window')
plt.savefig('fullmz_sn.pdf')


fig, ax = plt.subplots()
#markerline1, stemlines1, _ = ax.stem(df['m/z'],df['S/N'], linefmt='-',  markerfmt=" ",use_line_collection =True, label='150-2000 m/z')
#plt.setp(stemlines1, 'color', 'C0', 'linewidth', 2)
markerline2, stemlines2, _ = ax.stem(df2['m/z'],df2['S/N'], linefmt='-',  markerfmt=" ",use_line_collection =True, label='unique to 400-500 m/z')
plt.setp(stemlines2, 'color', 'C1', 'linewidth', 2)
ax.set_xlim(xlow, xhigh)
ax.set_ylim(0,ymax)
ax.set_xlabel('m/z')
ax.set_ylabel('S/N')
ax.set_title('400-500 m/z scan window')
plt.savefig('400-500mz_sn.pdf')

plt.close('all')



results_df = blank_subtracted_full
case = 'fullms_150-2000'

#results_df = narrow_unique
#case = 'unique_to_400-500mz'

filtered_results=results_df[(results_df['m/z']<800) & (results_df['S/N']>3)]

filtered_results['N']=filtered_results['N'].fillna(0)
filtered_results['O']=filtered_results['O'].fillna(0)
filtered_results['S']=filtered_results['S'].fillna(0)
filtered_results['P']=filtered_results['P'].fillna(0)
#filtered_results['Na']=filtered_results['Na'].fillna(0)

filtered_results['mol_class']='Unassigned'
filtered_results['mol_class'][filtered_results['C']>0]='CHO'
filtered_results['mol_class'][(filtered_results['C']>0) & (filtered_results['N']>0.5)]='CHON'
filtered_results['mol_class'][(filtered_results['C']>0) & (filtered_results['S']>0.5)]='CHOS'
filtered_results['mol_class'][(filtered_results['C']>0) & (filtered_results['P']>0.5)]='CHOP'
#filtered_results['mol_class'][(filtered_results['C']>0) & (filtered_results['Na']>0.5)]='CHONa'
filtered_results['mol_class'][(filtered_results['C']>0) & (filtered_results['S']>0.5) & (filtered_results['N']>0.5)]='CHONS'
filtered_results['mol_class'][(filtered_results['C']>0) & (filtered_results['P']>0.5) & (filtered_results['N']>0.5)]='CHONP'
#filtered_results['mol_class'][(filtered_results['C']>0) & (filtered_results['Na']>0.5) & (filtered_results['N']>0.5)]='CHONNa'
#filtered_results['mol_class'][(filtered_results['C']>0) & (filtered_results['P']>0.5) & (filtered_results['Na']>0.5) & (filtered_results['N']>0.5)]='CHONPNa'


results= filtered_results[filtered_results['Is Isotopologue']==0]
results['O/C']=results['O']/results['C']
results['H/C']=results['H']/results['C']
results['N/C']=results['N']/results['C']


## assignment error distribution
fig, ax = plt.subplots()
for mol_class in sorted(results['mol_class'].unique()):

    counts, bins = np.histogram(np.asarray(results[results['mol_class']==mol_class]['m/z Error (ppm)']),bins = 75)

    ax.plot(bins[:-1], counts, label = mol_class)

ax.set_xlim(-0.5,0.5)
ax.legend(frameon=False)
plt.savefig(case + '_assignment_error.pdf')

## van Krevelen plots 
## O/C, H/C
fig, ax = plt.subplots()
p1 =ax.scatter(x=results['O/C'],y=results['H/C'],s=results['Peak Height']/4000,c=results['Time'],cmap='viridis')
ax.set(xlabel='O/C',ylabel='H/C')
ax.legend(frameon=False)
fig.colorbar(p1,ax=[ax],label='Time')
plt.savefig(case+'_OC_HC.pdf')

## O/C, N/C
fig, ax = plt.subplots()
p1 = ax.scatter(x=results['O/C'],y=results['N/C'],s=results['Peak Height']/4000,c=results['Time'],cmap='viridis')
ax.set(xlabel='O/C',ylabel='N/C')
fig.colorbar(p1,ax=[ax],label='Time')
plt.savefig(case+'_OC_NC.pdf')




# normalized S/N space
xlow = 449
xhigh = 449.5
ymax = 200
df = blank_subtracted_full
df2 = blank_subtracted_n

fig, ax = plt.subplots()
markerline1, stemlines1, _ = ax.stem(df['m/z'],df['S/N'], linefmt='-',  markerfmt=" ",use_line_collection =True, label='150-2000 m/z')
plt.setp(stemlines1, 'color', 'C0', 'linewidth', 2)
markerline2, stemlines2, _ = ax.stem(df2['m/z'],df2['S/N'], linefmt='-',  markerfmt=" ",use_line_collection =True, label='unique to 400-500 m/z')
plt.setp(stemlines2, 'color', 'C1', 'linewidth', 2)
ax.set_xlim(xlow, xhigh)
ax.set_ylim(0,ymax)
ax.set_xlabel('m/z')
ax.set_ylabel('S/N')
ax.set_title('150-2000 m/z scan window')
plt.show()
plt.close('all')



sn_lim = 20
df2=narrow_unique[narrow_unique['S/N']<sn_lim]

df1=blank_subtracted_full[blank_subtracted_full['S/N']<sn_lim]
df2=narrow_unique[narrow_unique['S/N']<sn_lim]
df3=overlap[overlap['S/N']<sn_lim]

df2a = blank_subtracted_n
df1a = blank_subtracted_full

plot_ms(df2a,df2=df1a, start_mz=400,end_mz = 500, tstart = 10, lbls=['400-500 m/z', '150-2000 m/z']) 

plt.show()
plt.close('all')





def plot_ms(df1, start_mz, end_mz, tstart, df2=None,df3=None, assignment= None, ax_ms=None, lbls=None, norm=False, labs=False):   
    if ax_ms == None:
        f, ax = plt.subplots()
    
    else:
        ax = ax_ms

    mzrange= end_mz - start_mz
    ms_t_int=df1[df1['Time'] == tstart]
    ms_df = ms_t_int[((ms_t_int['Calibrated m/z']-start_mz)<mzrange) & ((ms_t_int['Calibrated m/z']-start_mz)>0)]

    maxdf1 = max(ms_df['S/N'])

    if norm:
        ms_df['S/N Norm'] = ms_df['S/N'] / maxdf1
    else:
        ms_df['S/N Norm'] = ms_df['S/N'] 

    

    print(ms_df['S/N Norm'])

    if lbls is not None:
        labels = lbls
    else:
        lbls = [None, None, None]

 
    _, stemlines1, _ =ax.stem('Calibrated m/z','S/N Norm',data=ms_df,  markerfmt=' ', basefmt=' ', linefmt='C0', label = labels[0])
    
    if df2 is not None:
        ms_t_int2=df2[df2['Time'] == tstart]
        ms_df2 = ms_t_int2[(abs(ms_t_int2['Calibrated m/z']-start_mz)<mzrange)& ((ms_t_int2['Calibrated m/z']-start_mz)>0)]

        maxdf2 = max(ms_df2['S/N'])

        if norm:
            ms_df2['S/N Norm'] = ms_df2['S/N'] / maxdf2
        else:
            ms_df2['S/N Norm'] = ms_df2['S/N'] 
        
        _, stemlines2, _ =ax.stem('Calibrated m/z','S/N Norm',data=ms_df2,  markerfmt=' ', basefmt=' ', linefmt='C1', label = labels[1])

    if df3 is not None:
        ms_t_int3=df3[df3['Time'] == tstart]
        ms_df3 = ms_t_int3[(abs(ms_t_int3['Calibrated m/z']-start_mz)<mzrange)& ((ms_t_int3['Calibrated m/z']-start_mz)>0)]

        maxdf3 = max(ms_df3['S/N'])

        if norm:
            ms_df3['S/N Norm'] = ms_df3['S/N'] / maxdf3
        else:
            ms_df3['S/N Norm'] = ms_df3['S/N'] 
        
        _, stemlines3, _ =ax.stem('Calibrated m/z','S/N Norm',data=ms_df3,  markerfmt=' ', basefmt=' ', linefmt='C2', label = labels[2])
    
    if df3 is not None:
        ax.set_ylim(0, max([maxdf1, maxdf2, maxdf3]) * 1.1)
    elif df2 is not None:
        ax.set_ylim(0, max([maxdf1, maxdf2]) * 1.1)
    else: 
        ax.set_ylim(0, maxdf1 * 1.1)

    ax.set_xlim(left = start_mz - mzrange*0.1, right = start_mz + mzrange + mzrange*0.1) 

    if labs:
        for mzr,peakr,mf,er in zip(ms_df['Calibrated m/z'], ms_df['S/N Norm'], ms_df['Molecular Formula'],  ms_df['m/z Error (ppm)']):

            #if (mzr- target_mz)  == 0:
            #    mz_text = ' m/z\n%.4f' % (mzr)
            #else:
            #    mz_text = r'$\Delta$' + ' m/z\n%.4f' % (mzr- target_mz)

            mz_text = 'm/z %.4f\n%s\n%.3f ppm' % (mzr,mf,er)
            ax.text(mzr, peakr + 0.02 *max(ms_df['S/N Norm']), mz_text, ha = 'center', fontsize = 'xx-small', weight = 'bold', color='C0')

        if df2 is not None:

            for mzr,peakr,mf, er in zip(ms_df2['Calibrated m/z'], ms_df2['S/N Norm'], ms_df2['Molecular Formula'], ms_df2['m/z Error (ppm)']):

            #if (mzr- target_mz)  == 0:
            #    mz_text = ' m/z\n%.4f' % (mzr)
            #else:
            #    mz_text = r'$\Delta$' + ' m/z\n%.4f' % (mzr- target_mz)

                mz_text = 'm/z %.4f\n%s\n%.3f ppm' % (mzr,mf,er)
                ax.text(mzr, peakr + 0.02 *max(ms_df2['S/N Norm']), mz_text, ha = 'center', fontsize = 'xx-small', weight = 'bold', color = 'C1')

        if df3 is not None:

            for mzr,peakr,mf, er in zip(ms_df3['Calibrated m/z'], ms_df3['S/N Norm'], ms_df3['Molecular Formula'], ms_df3['m/z Error (ppm)']):

            #if (mzr- target_mz)  == 0:
            #    mz_text = ' m/z\n%.4f' % (mzr)
            #else:
            #    mz_text = r'$\Delta$' + ' m/z\n%.4f' % (mzr- target_mz)

                mz_text = 'm/z %.4f\n%s\n%.3f ppm' % (mzr,mf,er)
                ax.text(mzr, peakr + 0.02 *max(ms_df3['S/N Norm']), mz_text, ha = 'center', fontsize = 'xx-small', weight = 'bold', color = 'C2')

   # theor_mz=pattern.mdiff+result['mass']
   # theor_int=pattern.ratio*result['abundance']
   # ax.stem(theor_mz,theor_int, basefmt=' ',linefmt='gray')

   # for isotope in pattern.isotope[pattern.requirement=='Y']:
   #     ax.stem('mz','intense',data=result[isotope],  markerfmt=' ', basefmt=' ',linefmt='red')
    if ax == None:

        ax.legend(('other', 'target'),bbox_to_anchor=(1.05, 1.0), loc='upper left',frameon=False)

        if(assignment):

            mf = assignment[0]
            score = assignment[1]
            er = assignment[2]

            ax.text(1.05,0.7,mf,transform=ax.transAxes)
            ax.text(1.05,0.6,'Error (ppm) = %.3f ' %er ,transform=ax.transAxes)
            ax.text(1.05,0.5,'Score = %.3f' %score ,transform=ax.transAxes)

    if norm: 
        ax.set(xlabel='Calibrated m/z',ylabel='Normalized S/N')
    else: 
        ax.set(xlabel='Calibrated m/z',ylabel='S/N')
    #ax.set_title('%.2f' %timerange[0] + ' to %.2f' %timerange[1] +' min', fontsize = 'medium')
    ax.legend(bbox_to_anchor = (1.05, 0.5), frameon =False, loc = 'center left')
    ax.axhline(y=0.0, color='black')
    plt.setp(stemlines1,'color', 'C0', 'linewidth', 2)
    if df2 is not None:
        plt.setp(stemlines2, 'color', 'C1','linewidth', 2)
    if df3 is not None:
        plt.setp(stemlines3, 'color', 'C2','linewidth', 2)
    plt.tight_layout()
    if ax_ms == None:
        return ax

