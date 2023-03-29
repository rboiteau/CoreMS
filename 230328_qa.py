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

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import seaborn as sns

import corems.lc_icpms_ftms.calc.lc_icrms_qc_assign as icrms
import corems.lc_icpms_ftms.calc.lc_icrms_helpers as lcmsfns
import corems.lc_icpms_ftms.calc.align_icpms_esi_fns as algnfns

"""
CoreMS run script for spring-env samples, collected at NHMFL in Nov 2023

Testing SN > 5, calorder = 1

Christian Dewey
28 Mar 23
"""

class Assignments:
    ''' 1 Mar 23
        Author: Christian Dewey
        
        intended for m/z windowing project
        works for 100 and 1000 m/z windows 
    '''
    def __init__(self, processed_assignments, mzrange, dataset_name):
        
        self._data_df = processed_assignments

        self.mzrange = mzrange

        self.dataset_name = dataset_name
        
        self.subsetByWindow()

        self.subsetAssignedUnique()

        self.combine()

        self.determineOverlap()


    def subsetByWindow(self):
        # create subsets of assignments based m/z window
        list100mz = []
        listfullzmz = [] 
        for mzwindow in self._data_df['Window Size (m/z)'].unique():

            r_df = self._data_df[self._data_df['Window Size (m/z)'] == mzwindow]

            if (mzwindow =='100'):
                list100mz.append(r_df) 

            elif (mzwindow =='1000'):
                listfullzmz.append(r_df) 

        self.mz100 = pd.concat(list100mz,ignore_index=True)
        self.mzfull = pd.concat(listfullzmz,ignore_index=True)

        self.mzfull_range = self.mzfull[(self.mzfull['m/z']>=self.mzrange[0])]
        self.mzfull_range = self.mzfull_range[(self.mzfull_range['m/z']<=self.mzrange[1])]

    def subsetAssignedUnique(self):
        # create subsets of assigned features and unique features 

        self.assigned_100mz = self.mz100[~self.mz100['Molecular Formula'].isnull()]
        self.assigned_fullmz = self.mzfull[~self.mzfull['Molecular Formula'].isnull()]  

        self.unique_100mz = lcmsfns.getUniqueFeatures(self.assigned_100mz)  ######
        self.unique_fullmz = lcmsfns.getUniqueFeatures(self.assigned_fullmz) #####

        self.assigned_fullmz_range = self.assigned_fullmz[(self.assigned_fullmz['m/z']>=self.mzrange[0])]
        self.assigned_fullmz_range = self.assigned_fullmz_range[(self.assigned_fullmz_range['m/z']<=self.mzrange[1])]
        self.unique_fullmz_range = lcmsfns.getUniqueFeatures(self.assigned_fullmz_range) #####

        print('\n100 m/z window (narrow):')
        print('%s features total\n%s assigned (%.1f%%)\n%s unique' %(np.shape(self.mz100)[0],  np.shape(self.assigned_100mz)[0], 
                                                                     (np.shape(self.assigned_100mz)[0] / np.shape(self.mz100)[0] * 100), 
                                                                     np.shape(self.unique_100mz)[0]))
        
        print('\n100-1100 m/z window, between %s-%s m/z:' %(self.mzrange[0],self.mzrange[1]))
        print('%s features total\n%s assigned (%.1f%%)\n%s unique' %(np.shape(self.mzfull_range)[0],np.shape(self.assigned_fullmz_range)[0],  
                                                                     np.shape(self.assigned_fullmz_range)[0] / np.shape(self.mzfull_range)[0] * 100, 
                                                                     len(self.unique_fullmz_range)))
        
        print('\n100-1100 m/z window (full):')
        print('%s features total\n%s assigned (%.1f%%)\n%s unique' %(np.shape(self.mzfull)[0],np.shape(self.assigned_fullmz)[0],  
                                                                     (np.shape(self.assigned_fullmz)[0] / np.shape(self.mzfull)[0] * 100),
                                                                     np.shape(self.unique_fullmz)[0]))

    def combine(self):
        # get combo unique 
        self.combo_assigned = pd.concat([self.assigned_100mz, self.assigned_fullmz_range], ignore_index=True)
        self.combo_unique =  pd.concat([self.unique_100mz, self.unique_fullmz_range], ignore_index=True)
        mzfull_range = self.mzfull[(self.mzfull['Calibrated m/z']>=self.mzrange[0]) & (self.mzfull['Calibrated m/z']<=self.mzrange[1])]
        self.combo_all = pd.concat([self.mz100, mzfull_range], ignore_index=True)

        print('\nOverall summary:\n%s features total' %len(self.combo_all))
        print('%s assigned features' %len(self.combo_assigned))
        print('%s unique features' %len(self.combo_unique))
   

    def determineOverlap(self):
        dd = self.combo_unique

        self.fullscan =  dd[dd['Window Size (m/z)'] == '1000']
        self.narrowscan = dd[dd['Window Size (m/z)'] == '100']

        print(len(self.fullscan),  len(self.narrowscan))

        self.narrowscan['mf_t'] = self.narrowscan['Molecular Formula'] + '-time_'+self.narrowscan['Time'].map(str)    
        self.fullscan['mf_t'] = self.fullscan['Molecular Formula'] + '-time_'+ self.fullscan['Time'].map(str)   
       
        shared_mf_t = list(set(list(self.narrowscan['mf_t'])) & set(list(self.fullscan['mf_t'])))
        onlyfull_mf_t = list(set(list(self.fullscan['mf_t'])) ^ set(shared_mf_t))
        onlynarrow_mf_t  = list(set(list(self.narrowscan['mf_t'])) ^ set(shared_mf_t))

        print('\n%s shared features' %len(shared_mf_t))
        print('%s features only in full window' %len(onlyfull_mf_t))
        print('%s features only in narrow window' %len(onlynarrow_mf_t))

        self.features_in_narrow_only = self.narrowscan[(self.narrowscan['mf_t'].isin(onlynarrow_mf_t))]
        self.features_in_full_only = self.fullscan[(self.fullscan['mf_t'].isin(onlyfull_mf_t))]

        dd2 = pd.concat([self.narrowscan,self.fullscan])
        self.features_in_both = dd2[(dd2['mf_t'].isin(shared_mf_t))]

        self.features_not_in_full = dd2[(~dd2['mf_t'].isin(shared_mf_t))]

        self.features_in_narrow_only['FeatureIn'] = '100'
        self.features_in_full_only['FeatureIn'] ='1000'
        self.features_in_both['FeatureIn'] = '100,1000'


def postAssignProcess():

    ##add 'Window Size (m/z)', 'm/z window', 'Rep', 'mol_class' columns

    global results
    global heter

    results = lcmsfns.add_mzwindow_col(results)                 # adds 'Window Size (m/z)' and 'm/z window' columns
    results = lcmsfns.addRepCol(results)                        # adds 'Rep' column
    molclasses = lcmsfns.get_mol_class(heter)                   # creates list of mol_classes based on heteroatom list
    results = lcmsfns.assign_mol_class(results,molclasses)      # adds 'mol_class' column 


def loadData():

    global data_dir
    global results 

    postAssignProcess()

    results.to_csv(data_dir+'p230328_spring-env_pos.csv')

    data_df = pd.read_csv(data_dir+'p230328_spring-env_pos.csv'  ) 

    data_df['Window Size (m/z)'] = data_df['Window Size (m/z)'].map(str)

    spring_asgn = Assignments(data_df, mzrange=[400,600], dataset_name='spring-env')

    return spring_asgn    


def plotCHO_CHON_mzerror(data_df):

    global plttitle, samplename, snl

    df_both = data_df.features_in_both
    df_full = pd.concat([data_df.features_in_full_only, df_both[df_both['Window Size (m/z)'] == '1000']])
    df_narrow = pd.concat([data_df.features_in_narrow_only, df_both[df_both['Window Size (m/z)'] == '100']])

    print('%s features in both windows' %len(df_both))
    print('%s features in full window' %len(df_full))
    print('%s features in narrow window' %len(df_narrow))

    molclass = ['CHON', 'CHO']

    xmin = -0.25
    xmax = 0.25

    hc = 'mol_class'

    fig, ((ax3,ax4,ax7),(ax5,ax6, ax8)) = plt.subplots(2, 3, figsize = (9,6))

    ax = ax3
    df = df_full[df_full['mol_class'].isin(molclass)]
    df = df[df['S/N']>snl]
    df = df.sort_values(by=['mol_class'])
    ax = sns.kdeplot(data=df,x='m/z Error (ppm)',hue=hc, ax=ax, legend = False, palette = sns.color_palette('colorblind'))
    ax.axvline(0, color = 'gray', linestyle = ':')
    ax.set_xlim(xmin,xmax)
    ax.set_title('%s features in full' %len(df), size = 8)

    ax = ax5
    df = df_narrow[df_narrow['mol_class'].isin(molclass)]
    df = df[df['S/N']>snl]
    df = df.sort_values(by=['mol_class'])
    ax = sns.kdeplot(data=df,x='m/z Error (ppm)',hue=hc, ax=ax, legend  = False, palette = sns.color_palette('colorblind'))
    ax.axvline(0, color = 'gray', linestyle = ':')
    ax.set_xlim(xmin,xmax)
    ax.set_xlim(xmin,xmax)
    ax.set_title('%s features in narrow' %len(df), size = 8)

    ax = ax4
    df = data_df.features_in_full_only[data_df.features_in_full_only['mol_class'].isin(molclass)]
    df = df[df['S/N']>snl]
    df = df.sort_values(by=['mol_class'])
    ax = sns.kdeplot(data=df,x='m/z Error (ppm)',hue=hc, ax=ax, legend = False, palette = sns.color_palette('colorblind'))
    ax.axvline(0, color = 'gray', linestyle = ':')
    ax.set_xlim(xmin,xmax)
    ax.set_title('%s features unique to full' %len(df), size = 8)

    ax = ax6
    df = data_df.features_in_narrow_only[data_df.features_in_narrow_only['mol_class'].isin(molclass)]
    df = df[df['S/N']>snl]
    df = df.sort_values(by=['mol_class'])
    ax = sns.kdeplot(data=df,x='m/z Error (ppm)',hue=hc, ax=ax, legend = False, palette = sns.color_palette('colorblind'))
    ax.axvline(0, color = 'gray', linestyle = ':')
    ax.set_xlim(xmin,xmax)
    ax.set_title('%s features unique to narrow' %len(df), size = 8)

    ax = ax7
    df = df_both[df_both['mol_class'].isin(molclass)]
    df = df[df['S/N']>snl]
    df = df.sort_values(by=['mol_class'])
    ax = sns.kdeplot(data=df,x='m/z Error (ppm)',hue=hc, ax=ax, palette = sns.color_palette('colorblind'))
    ax.axvline(0, color = 'gray', linestyle = ':')
    ax.set_xlim(xmin,xmax)
    ax.set_title('%s features in both' %len(df), size = 8)

    ax8.axis('off')

    lbls_art = []
    lbls = ['a','b','c','d','e']
    axs = [ax3,ax4,ax5,ax6, ax7]
    for lbl, ax in zip(lbls,axs):
        l = ax.text(-0.15, 1.05,lbl,
            horizontalalignment='center',
            verticalalignment='center',
            transform = ax.transAxes, fontweight='bold', fontsize = 10)
        lbls_art.append(l)

    legend = ax7.get_legend()
    handles = legend.legendHandles
    legend.remove()
    labels = sorted(molclass)
    leg = fig.legend(handles, labels, bbox_to_anchor=(0.8, 0.25), loc='center left',frameon=False, borderaxespad=0, title = 'Mol. Class', prop={'size': 8})

    sttl = fig.suptitle(plttitle)

    fig.tight_layout()

    sns.despine()
    lbla, lblb, lblc, lbld, lble= lbls_art
    flbl = plttitle.replace(" ","")
    if "/" in flbl:
        flbl = flbl.replace("/","")
    plt.savefig(data_dir + 'mz_error_CHO_CHON_' + flbl + '_' + samplename + '.pdf', bbox_extra_artists=(leg,lbla,lblb,lblc,lbld,lble,sttl), bbox_inches='tight')


def plotCHOS_CHOP_mzerror(data_df):

    global plttitle, samplename, snl

    df_both = data_df.features_in_both
    df_full = pd.concat([data_df.features_in_full_only, df_both[df_both['Window Size (m/z)'] == '1000']])
    df_narrow = pd.concat([data_df.features_in_narrow_only, df_both[df_both['Window Size (m/z)'] == '100']])

    print('%s features in both windows' %len(df_both))
    print('%s features in full window' %len(df_full))
    print('%s features in narrow window' %len(df_narrow))

    molclass = ['CHOS', 'CHOP']

    xmin = -0.25
    xmax = 0.25

    hc = 'mol_class'

    fig, ((ax3,ax4,ax7),(ax5,ax6, ax8)) = plt.subplots(2, 3, figsize = (9,6))

    ax = ax3
    df = df_full[df_full['mol_class'].isin(molclass)]
    df = df[df['S/N']>snl]
    df = df.sort_values(by=['mol_class'])
    ax = sns.kdeplot(data=df,x='m/z Error (ppm)',hue=hc, ax=ax, legend = False, palette = sns.color_palette('colorblind'))
    ax.axvline(0, color = 'gray', linestyle = ':')
    ax.set_xlim(xmin,xmax)
    ax.set_title('%s features in full' %len(df), size = 8)

    ax = ax5
    df = df_narrow[df_narrow['mol_class'].isin(molclass)]
    df = df[df['S/N']>snl]
    df = df.sort_values(by=['mol_class'])
    ax = sns.kdeplot(data=df,x='m/z Error (ppm)',hue=hc, ax=ax, legend  = False, palette = sns.color_palette('colorblind'))
    ax.axvline(0, color = 'gray', linestyle = ':')
    ax.set_xlim(xmin,xmax)
    ax.set_xlim(xmin,xmax)
    ax.set_title('%s features in narrow' %len(df), size = 8)

    ax = ax4
    df = data_df.features_in_full_only[data_df.features_in_full_only['mol_class'].isin(molclass)]
    df = df[df['S/N']>snl]
    df = df.sort_values(by=['mol_class'])
    ax = sns.kdeplot(data=df,x='m/z Error (ppm)',hue=hc, ax=ax, legend = False, palette = sns.color_palette('colorblind'))
    ax.axvline(0, color = 'gray', linestyle = ':')
    ax.set_xlim(xmin,xmax)
    ax.set_title('%s features unique to full' %len(df), size = 8)

    ax = ax6
    df = data_df.features_in_narrow_only[data_df.features_in_narrow_only['mol_class'].isin(molclass)]
    df = df[df['S/N']>snl]
    df = df.sort_values(by=['mol_class'])
    ax = sns.kdeplot(data=df,x='m/z Error (ppm)',hue=hc, ax=ax, legend = False, palette = sns.color_palette('colorblind'))
    ax.axvline(0, color = 'gray', linestyle = ':')
    ax.set_xlim(xmin,xmax)
    ax.set_title('%s features unique to narrow' %len(df), size = 8)

    ax = ax7
    df = df_both[df_both['mol_class'].isin(molclass)]
    df = df[df['S/N']>snl]
    df = df.sort_values(by=['mol_class'])
    ax = sns.kdeplot(data=df,x='m/z Error (ppm)',hue=hc, ax=ax, palette = sns.color_palette('colorblind'))
    ax.axvline(0, color = 'gray', linestyle = ':')
    ax.set_xlim(xmin,xmax)
    ax.set_title('%s features in both' %len(df), size = 8)

    ax8.axis('off')

    lbls_art = []
    lbls = ['a','b','c','d','e']
    axs = [ax3,ax4,ax5,ax6, ax7]
    for lbl, ax in zip(lbls,axs):
        l = ax.text(-0.15, 1.05,lbl,
            horizontalalignment='center',
            verticalalignment='center',
            transform = ax.transAxes, fontweight='bold', fontsize = 10)
        lbls_art.append(l)

    legend = ax7.get_legend()
    handles = legend.legendHandles
    legend.remove()
    labels = sorted(molclass)
    leg = fig.legend(handles, labels, bbox_to_anchor=(0.8, 0.25), loc='center left',frameon=False, borderaxespad=0, title = 'Mol. Class', prop={'size': 8})

    sttl = fig.suptitle(plttitle)

    fig.tight_layout()

    sns.despine()
    lbla, lblb, lblc, lbld, lble= lbls_art
    flbl = plttitle.replace(" ","")
    if "/" in flbl:
        flbl = flbl.replace("/","")
    plt.savefig(data_dir + 'mz_error_CHOS_CHOP_' + flbl + '_' + samplename + '.pdf', bbox_extra_artists=(leg,lbla,lblb,lblc,lbld,lble,sttl), bbox_inches='tight')


def plotCHONa_CHOK_mzerror(data_df):

    global plttitle, samplename, snl

    df_both = data_df.features_in_both
    df_full = pd.concat([data_df.features_in_full_only, df_both[df_both['Window Size (m/z)'] == '1000']])
    df_narrow = pd.concat([data_df.features_in_narrow_only, df_both[df_both['Window Size (m/z)'] == '100']])

    print('%s features in both windows' %len(df_both))
    print('%s features in full window' %len(df_full))
    print('%s features in narrow window' %len(df_narrow))

    elements = ['Na', 'K']
    ecols = ['C0', 'C1']
    xmin = -0.25
    xmax = 0.25

    elhue = {e: c for e, c in zip(elements, ecols)}

    fig, ((ax3,ax4,ax7),(ax5,ax6, ax8)) = plt.subplots(2, 3, figsize = (9,6))

    ax = ax3

    for e in elements:
        df = df_full[df_full[e]>0]
        df = df[df['S/N']>snl]
        df = df.sort_values(by=[e])
        ax = sns.kdeplot(data=df,x='m/z Error (ppm)',color = elhue[e], ax=ax, legend = False, palette = sns.color_palette('colorblind'))
    ax.axvline(0, color = 'gray', linestyle = ':')
    ax.set_xlim(xmin,xmax)
    ax.set_title('%s features in full' %len(df), size = 8)

    ax = ax5
    for e in elements:
        df = df_narrow[df_narrow[e]>0]
        df = df[df['S/N']>snl]
        df = df.sort_values(by=[e])
        ax = sns.kdeplot(data=df,x='m/z Error (ppm)',color = elhue[e],ax=ax, legend  = False, palette = sns.color_palette('colorblind'))
    ax.axvline(0, color = 'gray', linestyle = ':')
    ax.set_xlim(xmin,xmax)
    ax.set_xlim(xmin,xmax)
    ax.set_title('%s features in narrow' %len(df), size = 8)

    ax = ax4
    for e in elements:
        df = data_df.features_in_full_only[data_df.features_in_full_only[e]>0]
        df = df[df['S/N']>snl]
        df = df.sort_values(by=[e])
        ax = sns.kdeplot(data=df,x='m/z Error (ppm)', ax=ax, color = elhue[e], legend = False, palette = sns.color_palette('colorblind'))
    ax.axvline(0, color = 'gray', linestyle = ':')
    ax.set_xlim(xmin,xmax)
    ax.set_title('%s features unique to full' %len(df), size = 8)

    ax = ax6
    for e in elements:
        df = data_df.features_in_narrow_only[data_df.features_in_narrow_only[e]>0]
        df = df[df['S/N']>snl]
        df = df.sort_values(by=[e])
        ax = sns.kdeplot(data=df,x='m/z Error (ppm)',color = elhue[e], ax=ax, legend = False, palette = sns.color_palette('colorblind'))
    ax.axvline(0, color = 'gray', linestyle = ':')
    ax.set_xlim(xmin,xmax)
    ax.set_title('%s features unique to narrow' %len(df), size = 8)

    ax = ax7
    for e in elements:
        df = df_both[df_both[e]>0]
        df = df[df['S/N']>snl]
        df = df.sort_values(by=[e])
        ax = sns.kdeplot(data=df,x='m/z Error (ppm)',color = elhue[e], ax=ax, palette = sns.color_palette('colorblind'))
    ax.axvline(0, color = 'gray', linestyle = ':')
    ax.set_xlim(xmin,xmax)
    ax.set_title('%s features in both' %len(df), size = 8)

    ax8.axis('off')

    lbls_art = []
    lbls = ['a','b','c','d','e']
    axs = [ax3,ax4,ax5,ax6, ax7]
    for lbl, ax in zip(lbls,axs):
        l = ax.text(-0.15, 1.05,lbl,
            horizontalalignment='center',
            verticalalignment='center',
            transform = ax.transAxes, fontweight='bold', fontsize = 10)
        lbls_art.append(l)

    '''legend = ax7.get_legend()
    handles = legend.legendHandles
    legend.remove()
    labels = sorted(elements)
    leg = fig.legend(handles, labels, bbox_to_anchor=(0.8, 0.25), loc='center left',frameon=False, borderaxespad=0, title = 'Mol. Class', prop={'size': 8})
    '''
    sttl = fig.suptitle(plttitle)

    fig.tight_layout()

    sns.despine()
    lbla, lblb, lblc, lbld, lble= lbls_art
    flbl = plttitle.replace(" ","")
    if "/" in flbl:
        flbl = flbl.replace("/","")
    plt.savefig(data_dir + 'mz_error_CHONa_CHOK_' + flbl + '_' + samplename + '.pdf', bbox_extra_artists=(lbla,lblb,lblc,lbld,lble,sttl), bbox_inches='tight')


def plotCHOCu_CHOFe_mzerror(data_df):

    global plttitle, samplename, snl

    df_both = data_df.features_in_both
    df_full = pd.concat([data_df.features_in_full_only, df_both[df_both['Window Size (m/z)'] == '1000']])
    df_narrow = pd.concat([data_df.features_in_narrow_only, df_both[df_both['Window Size (m/z)'] == '100']])

    print('%s features in both windows' %len(df_both))
    print('%s features in full window' %len(df_full))
    print('%s features in narrow window' %len(df_narrow))

    elements = ['Fe']
    ecols = ['C0', 'C1']
    xmin = -0.25
    xmax = 0.25

    elhue = {e: c for e, c in zip(elements, ecols)}

    fig, ((ax3,ax4,ax7),(ax5,ax6, ax8)) = plt.subplots(2, 3, figsize = (9,6))

    ax = ax3

    for e in elements:
        df = df_full[df_full[e]>0]
        df = df[df['S/N']>snl]
        df = df.sort_values(by=[e])
        ax = sns.kdeplot(data=df,x='m/z Error (ppm)',color = elhue[e], ax=ax, legend = False, palette = sns.color_palette('colorblind'))
    ax.axvline(0, color = 'gray', linestyle = ':')
    ax.set_xlim(xmin,xmax)
    ax.set_title('%s features in full' %len(df), size = 8)

    ax = ax5
    for e in elements:
        df = df_narrow[df_narrow[e]>0]
        df = df[df['S/N']>snl]
        df = df.sort_values(by=[e])
        ax = sns.kdeplot(data=df,x='m/z Error (ppm)',color = elhue[e],ax=ax, legend  = False, palette = sns.color_palette('colorblind'))
    ax.axvline(0, color = 'gray', linestyle = ':')
    ax.set_xlim(xmin,xmax)
    ax.set_xlim(xmin,xmax)
    ax.set_title('%s features in narrow' %len(df), size = 8)

    ax = ax4
    for e in elements:
        df = data_df.features_in_full_only[data_df.features_in_full_only[e]>0]
        df = df[df['S/N']>snl]
        df = df.sort_values(by=[e])
        ax = sns.kdeplot(data=df,x='m/z Error (ppm)', ax=ax, color = elhue[e], legend = False, palette = sns.color_palette('colorblind'))
    ax.axvline(0, color = 'gray', linestyle = ':')
    ax.set_xlim(xmin,xmax)
    ax.set_title('%s features unique to full' %len(df), size = 8)

    ax = ax6
    for e in elements:
        df = data_df.features_in_narrow_only[data_df.features_in_narrow_only[e]>0]
        df = df[df['S/N']>snl]
        df = df.sort_values(by=[e])
        ax = sns.kdeplot(data=df,x='m/z Error (ppm)',color = elhue[e], ax=ax, legend = False, palette = sns.color_palette('colorblind'))
    ax.axvline(0, color = 'gray', linestyle = ':')
    ax.set_xlim(xmin,xmax)
    ax.set_title('%s features unique to narrow' %len(df), size = 8)

    ax = ax7
    for e in elements:
        df = df_both[df_both[e]>0]
        df = df[df['S/N']>snl]
        df = df.sort_values(by=[e])
        ax = sns.kdeplot(data=df,x='m/z Error (ppm)',color = elhue[e], ax=ax, palette = sns.color_palette('colorblind'))
    ax.axvline(0, color = 'gray', linestyle = ':')
    ax.set_xlim(xmin,xmax)
    ax.set_title('%s features in both' %len(df), size = 8)

    ax8.axis('off')

    lbls_art = []
    lbls = ['a','b','c','d','e']
    axs = [ax3,ax4,ax5,ax6, ax7]
    for lbl, ax in zip(lbls,axs):
        l = ax.text(-0.15, 1.05,lbl,
            horizontalalignment='center',
            verticalalignment='center',
            transform = ax.transAxes, fontweight='bold', fontsize = 10)
        lbls_art.append(l)

    '''legend = ax7.get_legend()
    handles = legend.legendHandles
    legend.remove()
    labels = sorted(elements)
    leg = fig.legend(handles, labels, bbox_to_anchor=(0.8, 0.25), loc='center left',frameon=False, borderaxespad=0, title = 'Mol. Class', prop={'size': 8})
    '''
    sttl = fig.suptitle(plttitle)

    fig.tight_layout()

    sns.despine()
    lbla, lblb, lblc, lbld, lble= lbls_art
    flbl = plttitle.replace(" ","")
    if "/" in flbl:
        flbl = flbl.replace("/","")
    plt.savefig(data_dir + 'mz_error_CHOCu_CHOFe_' + flbl + '_' + samplename + '.pdf', bbox_extra_artists=(lbla,lblb,lblc,lbld,lble,sttl), bbox_inches='tight')



if __name__ == '__main__':

    data_dir = '/home/dewey/Rawfiles/spring-env/pos/test/'
    samplename = 'spring-env-10-14min'
    fname = '230328_spring-env_pos.csv'
    heter = ['N', 'Na', 'S', 'P', 'K', 'Cu','Fe'] 
    results = pd.read_csv(data_dir+fname)
    plttitle = 'S/N > 3'
    snl = 3

    df = loadData()

    #plotCHO_CHON_mzerror(df)

    #plotCHOS_CHOP_mzerror(df)

    plotCHONa_CHOK_mzerror(df)

    plotCHOCu_CHOFe_mzerror(df)

    
