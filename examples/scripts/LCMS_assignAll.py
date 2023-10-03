import os
from tempfile import tempdir
from termios import TOSTOP
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

from corems.lc_icpms_ftms.factory.icpms_esims_helpers import extract, get_esiparser

import parMS


def plot_ms(masspectrum,srange,target_mz, assignment= None, ax_ms=None):   
    mass_spectrum = masspectrum
    if ax_ms == None:
        f, ax = plt.subplots()
    
    else:
        ax = ax_ms
    print('target mz: %s' %target_mz)

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

def subset_icpdata(element, icp_data, offset, timerange):
    etime = 'Time ' + element
    icp_subset = icp_data[[element,etime]]
    icp_subset[etime] = (icp_subset[etime] + offset) / 60 
    icp_subset = icp_subset[icp_subset[etime].between(timerange[0],timerange[1])]
    return icp_subset



def correlate(icpdata, EICdic, tic_times, offset, timerange):

    #EICdic = {}
    #pbar = tqdm.tqdm(mass_spectrum.mz_exp, desc="Getting EICs")
    
    #for mz in pbar:   
        #   print('AverageMS mz:' + str(mz))
    #    EIC=esi_parser.get_eics(target_mzs=[mz],tic_data={},peak_detection=False,smooth=False)
    #    EICdic[mz]=EIC[0][mz]
    timestart = timerange[0]
    timestop = timerange[1]
    times=tic_times
    
    icpsubset2 = subset_icpdata(heteroAtom,icpdata, offset = offset, timerange = [timestart, timestop])

    pbar = tqdm.tqdm(times, desc="Subsetting ICPMS data" )
    
    for i in pbar:
        icpsubset2.loc[-1]=['NaN',i]
        icpsubset2 = icpsubset2.sort_index().reset_index(drop=True)

    etime = 'Time ' + heteroAtom
    icpsubset2=icpsubset2.sort_values(by=etime)
    icpsubset2=icpsubset2.astype(float)
    icpsubset3=icpsubset2.interpolate()
    
    icp_interp=pd.DataFrame()
    pbar = tqdm.tqdm(times,desc="Interpolating ICPMS data")
    
    for i in pbar:
        icp_interp=icp_interp.append(icpsubset3[icpsubset3[etime]==i])

    mscorr={}
        #EICcurr=pd.DataFrame(index=icp_interp[etime],columns=['EIC',element])

    EICcurr=pd.DataFrame(index=icp_interp[etime],columns=['EIC',heteroAtom])
    EICcurr[heteroAtom]=icp_interp[heteroAtom].array
    
    pbar = tqdm.tqdm(EICdic.keys(),desc="Running correlation")
    
    for mz in pbar:
        EIC=pd.DataFrame({'EIC':EICdic[mz].eic,'Time':EICdic[mz].time})
        EIC_sub=EIC[EIC['Time'].between(timestart,timestop)]

        EICcurr['EIC']=EIC_sub['EIC'].values

        corvalue=EICcurr.corr(method='pearson')
        mscorr[mz]=corvalue.EIC[heteroAtom]**2   

        mzs_corr = pd.DataFrame.from_dict(mscorr,orient='index',columns=['corr'])

    return mzs_corr


class ScanList:
    
    def __init__(self, esi_parser, icpdata, chromaWindow, timerange):

        self.esi_parser = esi_parser 
        self.icpdata = icpdata 
        self.chroma_window = chromaWindow
        self.timerange = timerange

        self.getTICdf()
        self.getChromSub()
        self.getSrange()

    def getTICdf(self):
        self.tic = self.esi_parser.get_tic(ms_type = 'MS')[0]

        tic = self.tic

        scan_index = np.column_stack((range(len(tic.scans)), tic.scans, tic.time))

        self.scan_index = scan_index

    def getChromSub(self):

        tic = self.tic
        scan_freq  = tic.time[10] - tic.time[9] # min

        timestart = self.timerange[0]
        timestop = self.timerange[1]
        
        self.scans_per_chroma = int ( self.chroma_window  / (scan_freq * 60) ) + 1

        scan_index = self.scan_index
        self.sub = scan_index[np.where((scan_index[:,2]>timestart) & (scan_index[:,2]<=timestop))]

    def getSrange(self):

        temp = list(self.sub[:,1])
        scanlist = [int(i) for i in temp]

        self.srange = range(0,len(scanlist),self.scans_per_chroma)
        self.scanlist = scanlist

    def getMassSpectrumInScanRange(self, scanrange):
        mass_spectrum = self.esi_parser.get_average_mass_spectrum_by_scanlist(scanrange)
        self.mass_spectrum = mass_spectrum
        return self.mass_spectrum

    def corr(self, icpdata, offset):
        
        self.mzs_corr = correlate(icpdata,self.mass_spectrum, offset)


def assign(mass_spectrum,  elementDict, specific_valence=None):
    mass_spectrum.molecular_search_settings.error_method = 'None'
    mass_spectrum.molecular_search_settings.min_ppm_error = -10
    mass_spectrum.molecular_search_settings.max_ppm_error = 10

    mass_spectrum.molecular_search_settings.url_database = None
    mass_spectrum.molecular_search_settings.min_dbe = 0
    mass_spectrum.molecular_search_settings.max_dbe = 100

    elementList = elementDict.keys()
    for e in elementList:
    
        mass_spectrum.molecular_search_settings.usedAtoms[e] = elementDict[e]

    previous_atoms = list(mass_spectrum.molecular_search_settings.usedAtoms.keys())

    for e in previous_atoms:
        
        if e not in elementList:

            mass_spectrum.molecular_search_settings.usedAtoms.pop(e)

    if specific_valence:

        hetero = specific_valence[0]
        val = specific_valence[1]

        for e in elementList:
        
            if e == hetero:

                mass_spectrum.molecular_search_settings.used_atom_valences[e] = (val)

    mass_spectrum.molecular_search_settings.isProtonated = True
    mass_spectrum.molecular_search_settings.isRadical = False
    mass_spectrum.molecular_search_settings.isAdduct = False
    mass_spectrum.molecular_search_settings.ion_charge = ioncharge
    
    SearchMolecularFormulas(mass_spectrum, first_hit=False).run_worker_mass_spectrum()
    
    mass_spectrum.molecular_search_settings.score_method = "prob_score"
    mass_spectrum.molecular_search_settings.output_score_method = "prob_score"
    mass_spectrum.percentile_assigned(report_error=True)

    #mass_spectrum_by_classes = HeteroatomsClassification(mass_spectrum, choose_molecular_formula=True)

    assignments=mass_spectrum.to_dataframe()
    assignments=assignments.sort_values(by=['m/z'])
    
    return assignments



def extractEICs(scanlist_obj):
    
    scan_range_eics = {}

    srange = scanlist_obj.srange
    scan_index = scanlist_obj.scan_index
    scans_per_chroma = scanlist_obj.scans_per_chroma
    scanlist = scanlist_obj.scanlist
    
    esi_parser =   scanlist_obj.esi_parser

    range_labels = []
    rmax = int(scan_index[scan_index[:,1] == scanlist[-1]][0,0] - scan_index[scan_index[:,1] == scanlist[0]][0,0])

    for r in srange: 
            
            if r == rmax:
                
                break

            if (r + scans_per_chroma) > rmax:
                
                final_r = rmax
            else:
            
                final_r = r+scans_per_chroma
            
            range_label = '%s-%s' % (scanlist[r], scanlist[final_r-1])

            range_labels.append(range_label)

            print('\n...grabbing eics in scans ' + range_label)
            scanrange = scanlist[r:final_r]

              

            mass_spectrum = scanlist_obj.getMassSpectrumInScanRange(scanrange)
            
            #mzexp = mass_spectrum.mz_exp
            #ms = [ esi_parser, mzexp ]

#            pool = multiprocessing.Pool()

            #outputs_async = pool.map_async(extract,ms )
            #EICdic = outputs_async.get()
            #EICdic = extract(ms)
            #pool.close()
            #pool.join()
            EICdic = {}
            pbar = tqdm.tqdm(mass_spectrum.mz_exp, desc="Getting EICs")
            
            for mz in pbar:   
                #   print('AverageMS mz:' + str(mz))
                EIC=esi_parser.get_eics(target_mzs=[mz],tic_data={},peak_detection=False,smooth=False)
                EICdic[mz]=EIC[0][mz]

            scan_range_eics[range_label] = EICdic

    return scan_range_eics


def process(elementDict = None, scanlist_obj = None,  eics = None, assign_flag=True, returnLabels = False, specific_valence =None):
      
    srange = scanlist_obj.srange
    scan_index = scanlist_obj.scan_index
    scans_per_chroma = scanlist_obj.scans_per_chroma
    scanlist = scanlist_obj.scanlist
    icpdata = scanlist_obj.icpdata
    tic = scanlist_obj.esi_parser.get_tic(ms_type = 'MS')[0]
    tic_df = scanlist_obj.tic_df = pd.DataFrame({'time': tic.time,'scan': tic.scans})
    tic_times = tic_df[tic_df.time.between(timestart,timestop)].time.tolist()
    scan_range_data = {} 
    range_labels = []
    rmax = int(scan_index[scan_index[:,1] == scanlist[-1]][0,0] - scan_index[scan_index[:,1] == scanlist[0]][0,0])

    for r in srange: 
        
        if r == rmax:
            
            break

        if (r + scans_per_chroma) > rmax:
            
            final_r = rmax
        else:
           
            final_r = r+scans_per_chroma
        
        range_label = '%s-%s' % (scanlist[r], scanlist[final_r-1])

        range_labels.append(range_label)

        print('\n...processing scans ' + range_label)
        scanrange = scanlist[r:final_r]
    
        #mass_spectrum = esi_parser.get_average_mass_spectrum_by_scanlist(scanrange)

        mass_spectrum = scanlist_obj.getMassSpectrumInScanRange(scanrange)

        if eics == None:

            scan_range_eics = extractEICs()
        
        else:

            scan_range_eics = eics

        EICdic = scan_range_eics[range_label]

        mzs_corr = correlate(icpdata, EICdic, tic_times, offset, timerange = scanlist_obj.timerange)
        
        
        if assign_flag:

            if specific_valence:
                assignments = assign(mass_spectrum,elementDict,specific_valence)
            else:
                assignments = assign(mass_spectrum,elementDict)

            holder = np.zeros((len(assignments['m/z']),2))
            
            for mz,row in zip(assignments['m/z'], range(len(assignments['m/z']))):

                holder[row,1] = mzs_corr[mzs_corr.index == mz].iloc[0]['corr']
                holder[row,0] = mz

            pdholder = pd.DataFrame(holder, columns = ['m/z', 'corr'])

            assignments.insert(4,'corr',pdholder['corr'].values)
            assignments.insert(5,'mz2',pdholder['m/z'].values)


            scan_range_data[range_label] = (mass_spectrum, assignments)

        else:

            scan_range_data[range_label] = (mass_spectrum, mzs_corr)

    if returnLabels:

        return scan_range_data, range_labels

    else:

        return scan_range_data


def process_mp(elementDict = None, scanrange = None, scanlist_obj = None,  eics = None, assign_flag=True,  specific_valence =None):
          
    icpdata = scanlist_obj.icpdata

    tic = scanlist_obj.esi_parser.get_tic(ms_type = 'MS')[0]
    tic_df = scanlist_obj.tic_df = pd.DataFrame({'time': tic.time,'scan': tic.scans})
    tic_times = tic_df[tic_df.time.between(timestart,timestop)].time.tolist()

    scan_range_data = [] 

    mass_spectrum = scanlist_obj.getMassSpectrumInScanRange(scanrange)


    mzs_corr = correlate(icpdata, eic, tic_times, offset, timerange = [timestart,timestop]) # scanlist_obj.timerange)

    print(np.shape(mzs_corr))
    
    if assign_flag:

        if specific_valence:
            assignments = assign(mass_spectrum,elementDict,specific_valence)
        else:
            assignments = assign(mass_spectrum,elementDict)
        print(np.shape(assignments))
        print(assignments)
        holder = np.zeros((len(assignments['m/z']),2))
        print(np.shape(holder))
        print(holder)
        for mz,row in zip(assignments['m/z'], range(len(assignments['m/z']))):

            holder[row,1] = mzs_corr[mzs_corr.index == mz].iloc[0]['corr']
            holder[row,0] = mz

        pdholder = pd.DataFrame(holder, columns = ['m/z', 'corr'])

        assignments.insert(4,'corr',pdholder['corr'].values)
        #assignments.insert(5,'mz2',pdholder['m/z'].values)


        scan_range_data.append((mass_spectrum, assignments))

    else:

        scan_range_data.append((mass_spectrum, mzs_corr))

    return scan_range_data




def getResults(esi_parser, icpdata, scan_range_data, title, show, assigned, filter_corr, filter_hetero, plt_time_limits = None, icpms_max = None):
    
    range_labels = scan_range_data.keys()
    
    bestresults_list = []
    data_list = []

    tic = esi_parser.get_tic(ms_type = 'MS')[0]
    scan_index = np.column_stack((range(len(tic.scans)), tic.scans, tic.time))
 
    for rlabel in range_labels:

        temp = rlabel.split('-')
        scanrange = [int(temp[0]), int(temp[1])]
        data = scan_range_data[rlabel][1]
     
        range_start_min = scan_index[np.where(scan_index[:,1] == scanrange[0])][0,2]

        range_end_min = scan_index[np.where(scan_index[:,1] == scanrange[1])][0,2]

        timerange_for_scans =  [range_start_min, range_end_min]
        
        if filter_hetero: 

            match = re.findall(r'[A-Za-z]+|\d+', heteroAtom)
            match_element = match[1]

            bestresults=data[data[ match_element ] >= 1 ]
            
            inds = []
            starti = data.columns.get_loc(match_element)

            for i in range(starti+1,len(data.columns)):
               
                col = data.columns[i]
                
                temp = re.findall(r'[A-Za-z]+|\d+', col)
                
                match_i = temp[1]
                an_i = temp[0]

                if (match_element == match_i) and (len(col) == len(heteroAtom)) and (an_i.isdecimal()):
                    
                    inds.append(col)          
            
            if len(inds) > 1:
                
                for i in range(1, len(inds)):
                    
                    temp = data[data[ inds[i]] >= 1 ]
                    bestresults = pd.concat([bestresults, temp])
        
        else:
            
            bestresults = data

        
        if filter_corr:
            
            bestresults=bestresults[(bestresults['corr']>=threshold)]
            bestresults.sort_values(by='corr',ascending=False, inplace=True)
        
        
        if bestresults.empty:
                
                continue

        elif show:
            for r in range(np.shape(bestresults)[0]):
                res = r

                if assigned:
                    mz = bestresults.iloc[res]['m/z']
                    mf = bestresults.iloc[res]['Molecular Formula']
                    cs = bestresults.iloc[res]['Confidence Score']
                    mz_er = bestresults.iloc[res]['m/z Error (ppm)']
                    assignment_list = [mf,cs,mz_er]

                else:
                    mz = bestresults.index.values[res]

                cor_mz = bestresults.iloc[res]['corr']
                
                eicb=esi_parser.get_eics(target_mzs=[mz],tic_data={},peak_detection=False,smooth=False)
                eict=eicb[0][mz].time
                eics=eicb[0][mz].eic
                icpt = (icpdata['Time ' + heteroAtom] + offset) /60                                            
                icps = icpdata[heteroAtom]

            
                fig, ax = plt.subplots()

                ax2 = ax.twinx()

                ax.plot(icpt,icps, 'g-')

                ax2.plot(eict, eics, 'b-')

                ax.set_xlabel('Time (min)')
                ax.set_ylabel('%s ICP-MS signal' % heteroAtom, color='g')
                ax2.set_ylabel('ESI-MS signal', color='b')

                ax.yaxis.label.set_color('g')
                ax.tick_params(axis='y', colors='g')

                ax2.yaxis.label.set_color('b')
                ax2.tick_params(axis='y', colors='b')
                
                if plt_time_limits and (icpms_max == None):
                    ax.set_xlim(plt_time_limits[0], plt_time_limits[1])
                    icpsub_lim = subset_icpdata(icp_data = icpdata, offset=0,timerange=[plt_time_limits[0],plt_time_limits[1]],element=heteroAtom)
                    ymax =1.1*(max(icpsub_lim[heteroAtom]))
                    ax.set_ylim(-0.005*ymax, ymax)
                    ax2.set_ylim(-0.005*max(eics),max(eics)*1.1)

                elif plt_time_limits and (icpms_max != None):
                    ax.set_xlim(plt_time_limits[0], plt_time_limits[1])
                    ax.set_ylim(-0.005*icpms_max, icpms_max)
                    ax2.set_ylim(-0.005*max(eics),max(eics)*1.1)

                
                elif plt_time_limits == None:
                    ax.set_xlim([timestart, timestop])

                ax.set_title(title)

                ax.text(0.05, 0.8, '$R^2$: %.3f' % cor_mz,transform=ax.transAxes)

                if assigned and (pd.isna(mf) is False):    
                    ax.text(0.05, 0.95, 'm/z: %.4f' %mz,transform=ax.transAxes)
                    ax.text(0.05, 0.90, 'Formula: %s' %(mf),transform=ax.transAxes)
                    ax.text(0.05,0.75, 'm/z error (ppm): %.3f' %mz_er,transform=ax.transAxes)
                    ax.text(0.05, 0.675, 'Confidence: %.3f' %cs,transform=ax.transAxes)

                elif assigned and (pd.isna(mf) is True):    
                    ax.text(0.05, 0.95, 'm/z: %.4f' %mz,transform=ax.transAxes)
                    ax.text(0.05, 0.875, 'Formula not assigned',transform=ax.transAxes)

                plt.show()
                ms_i = scan_range_data[rlabel][0]
                ax = plot_ms(masspectrum = ms_i, srange=[scanrange, timerange_for_scans], target_mz=mz, assignment=assignment_list)
                
                plt.show()
           
        else:

            nplts = np.shape(bestresults)[0]
            
            if nplts > 1:
                
                fig = plt.figure(figsize=(8, 4*nplts ))
                gs = fig.add_gridspec(ncols=2, nrows=nplts, width_ratios=[1.2, 1],hspace = hspace1)
                fig.suptitle(title,va = 'center')
                        
                
                for r in range(0,nplts):
                #print(axs)
                #for r, ax in zip(range(nplts), axs):
                    ax1 = fig.add_subplot(gs[2*r])
                    
                    #print(axs)
                    res = r

                    if assigned:
                        
                        mz = bestresults.iloc[res]['m/z']
                        mf = bestresults.iloc[res]['Molecular Formula']
                        cs = bestresults.iloc[res]['Confidence Score']
                        mz_er = bestresults.iloc[res]['m/z Error (ppm)']
                        assignment_list = [mf,cs,mz_er]

                    else:
                        
                        mz = bestresults.index.values[res]

                    cor_mz = bestresults.iloc[res]['corr']
                    
                    eicb=esi_parser.get_eics(target_mzs=[mz],tic_data={},peak_detection=False,smooth=False)
                    eict=eicb[0][mz].time
                    eics=eicb[0][mz].eic
                    icpt = (icpdata['Time ' + heteroAtom] + offset) /60
                    icps = icpdata[heteroAtom]

                    #ax1 = ax
                    ax2 = ax1.twinx()

                    line1, = ax1.plot(icpt,icps, 'g-', linewidth = 1, label = heteroAtom)
                    line2, = ax2.plot(eict, eics, 'b-.', linewidth = 0.8, label = 'EIC %.4f' %mz)

                    ax1.set_xlabel('Time (min)')
                    ax1.set_ylabel('%s ICP-MS signal intensity' % heteroAtom, color='g')
                    ax2.set_ylabel('ESI-MS signal intensity', color='b')

                    ax1.yaxis.label.set_color('g')
                    ax1.tick_params(axis='y', colors='g')

                    ax2.yaxis.label.set_color('b')
                    ax2.tick_params(axis='y', colors='b')

                    
                    if plt_time_limits and (icpms_max == None):
                        ax1.set_xlim(plt_time_limits[0], plt_time_limits[1])
                        icpsub_lim = subset_icpdata(icp_data = icpdata, offset=0,timerange=[plt_time_limits[0],plt_time_limits[1]],element=heteroAtom)
                        ymax =1.1*(max(icpsub_lim[heteroAtom]))
                        ax1.set_ylim(-0.005*ymax, ymax)
                        ax2.set_ylim(-0.005*max(eics),max(eics)*1.1)

                    elif plt_time_limits and (icpms_max != None):
                        ax1.set_xlim(plt_time_limits[0], plt_time_limits[1])
                        ax1.set_ylim(-0.005*icpms_max, icpms_max)
                        ax2.set_ylim(-0.005*max(eics),max(eics)*1.1)
                    
                    elif plt_time_limits == None:
                        ax1.set_xlim([timestart, timestop])

                    ax1.text(0.05, 0.90, '$R^2$: %.3f' % cor_mz,transform=ax1.transAxes, fontsize = 'small')

                    if assigned and (pd.isna(mf) is False):    
                        ax1.text(0.05, 0.95, 'm/z: %.6f' %mz,transform=ax1.transAxes, fontsize = 'small')
                        ax1.text(0.05, 0.85, 'Formula: %s' %(mf),transform=ax1.transAxes, fontsize = 'small')
                        ax1.text(0.05,0.8, 'm/z error (ppm): %.3f' %mz_er,transform=ax1.transAxes, fontsize = 'small')
                        ax1.text(0.05, 0.75, 'Confidence: %.3f' %cs,transform=ax1.transAxes, fontsize = 'small')

                    elif assigned and (pd.isna(mf) is True):    
                        ax1.text(0.05, 0.95, 'm/z: %.6f' %mz,transform=ax1.transAxes, fontsize = 'small')
                        ax1.text(0.05, 0.85, 'Formula not assigned',transform=ax1.transAxes, fontsize = 'small')

                    ax1.legend(handles=[line1, line2], loc = legendloc, frameon = False, fontsize = 'small')   

                    ms_i = scan_range_data[rlabel][0]
                
                    ax_ms = fig.add_subplot(gs[2*r + 1])  
                    plot_ms(masspectrum = ms_i, srange=[scanrange, timerange_for_scans], target_mz=mz, assignment=assignment_list, ax_ms=ax_ms)

                    plt.subplots_adjust(bottom=bplt1, left = lplt1, right=rplt1, top=tplt1, hspace=hspace1, wspace =wspace1)
                
                # plt.subplots_adjust(bottom=bplt1, left = lplt1, right=rplt1, top=tplt1, hspace=hspace1, wspace =wspace1)

            else:

                fig, axb = plt.subplots(figsize=(8, 4 ))

                axb.set_title(title)
                
                res = 0

                if assigned:
                    
                    mz = bestresults.iloc[res]['m/z']
                    mf = bestresults.iloc[res]['Molecular Formula']
                    cs = bestresults.iloc[res]['Confidence Score']
                    mz_er = bestresults.iloc[res]['m/z Error (ppm)']
                    assignment_list = [mf,cs,mz_er]

                else:
                    
                    mz = bestresults.index.values[res]

                cor_mz = bestresults.iloc[res]['corr']
                
                eicb=esi_parser.get_eics(target_mzs=[mz],tic_data={},peak_detection=False,smooth=False)
                eict=eicb[0][mz].time
                eics=eicb[0][mz].eic
                icpt = (icpdata['Time ' + heteroAtom] + offset) /60
                icps = icpdata[heteroAtom]

                ax1 = axb
                ax2 = ax1.twinx()

                ax1.plot(icpt,icps, 'g-', linewidth = 1)
                ax2.plot(eict, eics, 'b-.', linewidth = 0.8)

                ax1.set_xlabel('Time (min)')
                ax1.set_ylabel('%s ICP-MS signal' % heteroAtom, color='g')
                ax2.set_ylabel('ESI-MS signal', color='b')

                ax1.yaxis.label.set_color('g')
                ax1.tick_params(axis='y', colors='g')

                ax2.yaxis.label.set_color('b')
                ax2.tick_params(axis='y', colors='b')
                
                if plt_time_limits and (icpms_max == None):
                    axb.set_xlim(plt_time_limits[0], plt_time_limits[1])
                    icpsub_lim = subset_icpdata(icp_data = icpdata, offset=0,timerange=[plt_time_limits[0],plt_time_limits[1]],element=heteroAtom)
                    ymax =1.1*(max(icpsub_lim[heteroAtom]))
                    axb.set_ylim(-0.005*ymax, ymax)
                    ax2.set_ylim(-0.005*max(eics),max(eics)*1.1)

                elif plt_time_limits and (icpms_max != None):
                    axb.set_ylim(-0.005*icpms_max, icpms_max)
                    ax2.set_ylim(-0.005*max(eics),max(eics)*1.1)
                    axb.set_xlim(plt_time_limits[0], plt_time_limits[1])
               
                elif plt_time_limits == None:
                    axb.set_xlim([timestart, timestop])

                axb.text(0.05, 0.80, '$R^2$: %.3f' % cor_mz,transform=axb.transAxes)
                
                if assigned and (pd.isna(mf) is False):    
                    axb.text(0.05, 0.95, 'm/z: %.6f' %mz,transform=axb.transAxes)
                    axb.text(0.05, 0.90, 'Formula: %s' %(mf),transform=axb.transAxes)
                    axb.text(0.05,0.75, 'm/z error (ppm): %.3f' %mz_er,transform=axb.transAxes)
                    axb.text(0.05, 0.675, 'Confidence: %.3f' %cs,transform=axb.transAxes)

                if assigned and (pd.isna(mf) is True):    
                    axb.text(0.05, 0.95, 'm/z: %.6f' %mz,transform=axb.transAxes)
                    axb.text(0.05, 0.875, 'Formula not assigned',transform=axb.transAxes)
            
                ms_i = scan_range_data[rlabel][0]
                ax = plot_ms(masspectrum = ms_i, srange=[scanrange, timerange_for_scans], target_mz=mz, assignment=assignment_list)

            dir = '/Users/christiandewey/Downloads/'
            
            plt.savefig(dir + title + '_' + rlabel + '.pdf',bbox_inches='tight')

            plt.close('all')

        bestresults_list.append(bestresults)
        data_list.append(data)

    return bestresults_list, data_list



def getResults_mp(esi_parser, icpdata, scan_range_data, title, show, assigned, filter_corr, filter_hetero, plt_time_limits = None, icpms_max = None):
    
    #range_labels = scan_range_data.keys()
    
    bestresults_list = []
    data_list = []
 
    #for rlabel in range_labels:
    #temp = rlabel.split('-')
    #scanrange = [int(temp[0]), int(temp[1])]
    data = scan_range_data[0][1]  # index 0 is tuple of ms and assignments 
    
    if filter_hetero: 

        match = re.findall(r'[A-Za-z]+|\d+', heteroAtom)
        match_element = match[1]

        bestresults=data[data[ match_element ] >= 1 ]
        
        inds = []
        starti = data.columns.get_loc(match_element)

        for i in range(starti+1,len(data.columns)):
            
            col = data.columns[i]
            
            temp = re.findall(r'[A-Za-z]+|\d+', col)
            
            match_i = temp[1]
            an_i = temp[0]

            if (match_element == match_i) and (len(col) == len(heteroAtom)) and (an_i.isdecimal()):
                
                inds.append(col)          
        
        if len(inds) > 1:
            
            for i in range(1, len(inds)):
                
                temp = data[data[ inds[i]] >= 1 ]
                bestresults = pd.concat([bestresults, temp])
    
    else:
        
        bestresults = data

    
    if filter_corr:
        
        bestresults=bestresults[(bestresults['corr']>=threshold)]
        bestresults.sort_values(by='corr',ascending=False, inplace=True)
    

    if show:
        for r in range(np.shape(bestresults)[0]):
            res = r

            if assigned:
                mz = bestresults.iloc[res]['m/z']
                mf = bestresults.iloc[res]['Molecular Formula']
                cs = bestresults.iloc[res]['Confidence Score']
                mz_er = bestresults.iloc[res]['m/z Error (ppm)']
                assignment_list = [mf,cs,mz_er]

            else:
                mz = bestresults.index.values[res]

            cor_mz = bestresults.iloc[res]['corr']
            
            eicb=esi_parser.get_eics(target_mzs=[mz],tic_data={},peak_detection=False,smooth=False)
            eict=eicb[0][mz].time
            eics=eicb[0][mz].eic
            icpt = (icpdata['Time ' + heteroAtom] + offset) /60                                            
            icps = icpdata[heteroAtom]

        
            fig, ax = plt.subplots()

            ax2 = ax.twinx()

            ax.plot(icpt,icps, 'g-')

            ax2.plot(eict, eics, 'b-')

            ax.set_xlabel('Time (min)')
            ax.set_ylabel('%s ICP-MS signal' % heteroAtom, color='g')
            ax2.set_ylabel('ESI-MS signal', color='b')

            ax.yaxis.label.set_color('g')
            ax.tick_params(axis='y', colors='g')

            ax2.yaxis.label.set_color('b')
            ax2.tick_params(axis='y', colors='b')
            
            if plt_time_limits and (icpms_max == None):
                ax.set_xlim(plt_time_limits[0], plt_time_limits[1])
                icpsub_lim = subset_icpdata(icp_data = icpdata, offset=0,timerange=[plt_time_limits[0],plt_time_limits[1]],element=heteroAtom)
                ymax =1.1*(max(icpsub_lim[heteroAtom]))
                ax.set_ylim(-0.005*ymax, ymax)
                ax2.set_ylim(-0.005*max(eics),max(eics)*1.1)

            elif plt_time_limits and (icpms_max != None):
                ax.set_xlim(plt_time_limits[0], plt_time_limits[1])
                ax.set_ylim(-0.005*icpms_max, icpms_max)
                ax2.set_ylim(-0.005*max(eics),max(eics)*1.1)

            
            elif plt_time_limits == None:
                ax.set_xlim([timestart, timestop])

            ax.set_title(title)

            ax.text(0.05, 0.8, '$R^2$: %.3f' % cor_mz,transform=ax.transAxes)

            if assigned and (pd.isna(mf) is False):    
                ax.text(0.05, 0.95, 'm/z: %.4f' %mz,transform=ax.transAxes)
                ax.text(0.05, 0.90, 'Formula: %s' %(mf),transform=ax.transAxes)
                ax.text(0.05,0.75, 'm/z error (ppm): %.3f' %mz_er,transform=ax.transAxes)
                ax.text(0.05, 0.675, 'Confidence: %.3f' %cs,transform=ax.transAxes)

            elif assigned and (pd.isna(mf) is True):    
                ax.text(0.05, 0.95, 'm/z: %.4f' %mz,transform=ax.transAxes)
                ax.text(0.05, 0.875, 'Formula not assigned',transform=ax.transAxes)

            plt.show()
            ms_i = scan_range_data[0][0]
            ax = plot_ms(masspectrum = ms_i, srange=[scanrange, [timestart,timestop]], target_mz=mz, assignment=assignment_list)
            
            plt.show()
        
    else:

        nplts = np.shape(bestresults)[0]
        print(nplts)

        if nplts > 1:
            print('if')
            fig = plt.figure(figsize=(8, 4*nplts ))
            gs = fig.add_gridspec(ncols=2, nrows=nplts, width_ratios=[1.2, 1],hspace = hspace1)
            fig.suptitle(title,va = 'center')
                    
            print(bestresults)
            for r in range(0,nplts):
            #print(axs)
            #for r, ax in zip(range(nplts), axs):
                ax1 = fig.add_subplot(gs[2*r])
                print(r)
                #print(axs)
                res = r

                if assigned:
                    
                    mz = bestresults.iloc[res]['m/z']
                    mf = bestresults.iloc[res]['Molecular Formula']
                    cs = bestresults.iloc[res]['Confidence Score']
                    mz_er = bestresults.iloc[res]['m/z Error (ppm)']
                    assignment_list = [mf,cs,mz_er]

                else:
                    
                    mz = bestresults.index.values[res]

                cor_mz = bestresults.iloc[res]['corr']
                
                eicb=esi_parser.get_eics(target_mzs=[mz],tic_data={},peak_detection=False,smooth=False)
                eict=eicb[0][mz].time
                eics=eicb[0][mz].eic
                icpt = (icpdata['Time ' + heteroAtom] + offset) /60
                icps = icpdata[heteroAtom]

                #ax1 = ax
                ax2 = ax1.twinx()

                line1, = ax1.plot(icpt,icps, 'g-', linewidth = 1, label = heteroAtom)
                line2, = ax2.plot(eict, eics, 'b-.', linewidth = 0.8, label = 'EIC %.4f' %mz)

                ax1.set_xlabel('Time (min)')
                ax1.set_ylabel('%s ICP-MS signal intensity' % heteroAtom, color='g')
                ax2.set_ylabel('ESI-MS signal intensity', color='b')

                ax1.yaxis.label.set_color('g')
                ax1.tick_params(axis='y', colors='g')

                ax2.yaxis.label.set_color('b')
                ax2.tick_params(axis='y', colors='b')

                
                if plt_time_limits and (icpms_max == None):
                    ax1.set_xlim(plt_time_limits[0], plt_time_limits[1])
                    icpsub_lim = subset_icpdata(icp_data = icpdata, offset=0,timerange=[plt_time_limits[0],plt_time_limits[1]],element=heteroAtom)
                    ymax =1.1*(max(icpsub_lim[heteroAtom]))
                    ax1.set_ylim(-0.005*ymax, ymax)
                    ax2.set_ylim(-0.005*max(eics),max(eics)*1.1)

                elif plt_time_limits and (icpms_max != None):
                    ax1.set_xlim(plt_time_limits[0], plt_time_limits[1])
                    ax1.set_ylim(-0.005*icpms_max, icpms_max)
                    ax2.set_ylim(-0.005*max(eics),max(eics)*1.1)
                
                elif plt_time_limits == None:
                    ax1.set_xlim([timestart, timestop])

                ax1.text(0.05, 0.90, '$R^2$: %.3f' % cor_mz,transform=ax1.transAxes, fontsize = 'small')

                if assigned and (pd.isna(mf) is False):    
                    ax1.text(0.05, 0.95, 'm/z: %.6f' %mz,transform=ax1.transAxes, fontsize = 'small')
                    ax1.text(0.05, 0.85, 'Formula: %s' %(mf),transform=ax1.transAxes, fontsize = 'small')
                    ax1.text(0.05,0.8, 'm/z error (ppm): %.3f' %mz_er,transform=ax1.transAxes, fontsize = 'small')
                    ax1.text(0.05, 0.75, 'Confidence: %.3f' %cs,transform=ax1.transAxes, fontsize = 'small')

                elif assigned and (pd.isna(mf) is True):    
                    ax1.text(0.05, 0.95, 'm/z: %.6f' %mz,transform=ax1.transAxes, fontsize = 'small')
                    ax1.text(0.05, 0.85, 'Formula not assigned',transform=ax1.transAxes, fontsize = 'small')

                ax1.legend(handles=[line1, line2], loc = legendloc, frameon = False, fontsize = 'small')   

                ms_i = scan_range_data[0][0]
            
                ax_ms = fig.add_subplot(gs[2*r + 1])  
                plot_ms(masspectrum = ms_i, srange=[scanrange, [timestart,timestop]], target_mz=mz, assignment=assignment_list, ax_ms=ax_ms)

                plt.subplots_adjust(bottom=bplt1, left = lplt1, right=rplt1, top=tplt1, hspace=hspace1, wspace =wspace1)
            
            # plt.subplots_adjust(bottom=bplt1, left = lplt1, right=rplt1, top=tplt1, hspace=hspace1, wspace =wspace1)

        else:

            fig, axb = plt.subplots(figsize=(8, 4 ))

            axb.set_title(title)
            
            res = 0

            if assigned:
                
                mz = bestresults.iloc[res]['m/z']
                mf = bestresults.iloc[res]['Molecular Formula']
                cs = bestresults.iloc[res]['Confidence Score']
                mz_er = bestresults.iloc[res]['m/z Error (ppm)']
                assignment_list = [mf,cs,mz_er]

            else:
                
                mz = bestresults.index.values[res]

            cor_mz = bestresults.iloc[res]['corr']
            
            eicb=esi_parser.get_eics(target_mzs=[mz],tic_data={},peak_detection=False,smooth=False)
            eict=eicb[0][mz].time
            eics=eicb[0][mz].eic
            icpt = (icpdata['Time ' + heteroAtom] + offset) /60
            icps = icpdata[heteroAtom]

            ax1 = axb
            ax2 = ax1.twinx()

            ax1.plot(icpt,icps, 'g-', linewidth = 1)
            ax2.plot(eict, eics, 'b-.', linewidth = 0.8)

            ax1.set_xlabel('Time (min)')
            ax1.set_ylabel('%s ICP-MS signal' % heteroAtom, color='g')
            ax2.set_ylabel('ESI-MS signal', color='b')

            ax1.yaxis.label.set_color('g')
            ax1.tick_params(axis='y', colors='g')

            ax2.yaxis.label.set_color('b')
            ax2.tick_params(axis='y', colors='b')
            
            if plt_time_limits and (icpms_max == None):
                axb.set_xlim(plt_time_limits[0], plt_time_limits[1])
                icpsub_lim = subset_icpdata(icp_data = icpdata, offset=0,timerange=[plt_time_limits[0],plt_time_limits[1]],element=heteroAtom)
                ymax =1.1*(max(icpsub_lim[heteroAtom]))
                axb.set_ylim(-0.005*ymax, ymax)
                ax2.set_ylim(-0.005*max(eics),max(eics)*1.1)

            elif plt_time_limits and (icpms_max != None):
                axb.set_ylim(-0.005*icpms_max, icpms_max)
                ax2.set_ylim(-0.005*max(eics),max(eics)*1.1)
                axb.set_xlim(plt_time_limits[0], plt_time_limits[1])
            
            elif plt_time_limits == None:
                axb.set_xlim([timestart, timestop])

            axb.text(0.05, 0.80, '$R^2$: %.3f' % cor_mz,transform=axb.transAxes)
            
            if assigned and (pd.isna(mf) is False):    
                axb.text(0.05, 0.95, 'm/z: %.6f' %mz,transform=axb.transAxes)
                axb.text(0.05, 0.90, 'Formula: %s' %(mf),transform=axb.transAxes)
                axb.text(0.05,0.75, 'm/z error (ppm): %.3f' %mz_er,transform=axb.transAxes)
                axb.text(0.05, 0.675, 'Confidence: %.3f' %cs,transform=axb.transAxes)

            if assigned and (pd.isna(mf) is True):    
                axb.text(0.05, 0.95, 'm/z: %.6f' %mz,transform=axb.transAxes)
                axb.text(0.05, 0.875, 'Formula not assigned',transform=axb.transAxes)
        
            ms_i = scan_range_data[0][0]
            ax = plot_ms(masspectrum = ms_i, srange=[scanrange, [timestart,timestop]], target_mz=mz, assignment=assignment_list)

        dir = '/Users/christiandewey/Downloads/'
        
        plt.savefig(dir + title + '.pdf',bbox_inches='tight')

        plt.close('all')

    #bestresults_list.append(bestresults)
    #data_list.append(data)

    return bestresults_list, data_list


def plot_single_mz(esi_parser,icpdata,scan_range_data_i, scanrange, mz, results = None, show = True, plt_time_limits = None, title = None, assigned = True, icpms_max = None):

    eicb=esi_parser.get_eics(target_mzs=[mz],tic_data={},peak_detection=False,smooth=False)
    print(scanrange)
    print(mz)

    if results is not None:
        
        for mzr in results['m/z']:
            
            mzr_str = '%.6f' % mzr
            print(mzr_str)
            if (mzr_str == str(mz)) | (mzr == mz):
                print('here')
                sub = results.loc[results['m/z'] == mzr]

        if assigned:
            mflist = []
            cslist = []
            erlist = []

            for i in range(len(sub.index)):
                mf = sub['Molecular Formula'].iloc[i]
                cs = sub['Confidence Score'].iloc[i]
                mz_er = sub['m/z Error (ppm)'].iloc[i]
                assignment_list = [mf,cs,mz_er]

                mflist.append(mf)
                cslist.append(cs)
                erlist.append(mz_er)

        cor_mz = sub['corr'].iloc[0]

    print('mf: %s' % mf)
    eict=eicb[0][mz].time
    eics=eicb[0][mz].eic
    icpt = (icpdata['Time ' + heteroAtom] + offset) /60
    icps = icpdata[heteroAtom]


    fig = plt.figure(figsize=(8, 4 ))
    gs = fig.add_gridspec(ncols=2, nrows=1, width_ratios=[1.2, 1],hspace = hspace1)
    fig.suptitle(title,va = 'center')


    ax = fig.add_subplot(gs[0])


    ax2 = ax.twinx()



    line1, = ax.plot(icpt,icps, 'g-', linewidth = 1, label = heteroAtom)
    line2, = ax2.plot(eict, eics, 'b-.', linewidth = 0.8, label = 'EIC %.4f' %mz)

    ax.set_xlabel('Time (min)')
    ax.set_ylabel('%s ICP-MS signal intensity' % heteroAtom, color='g')
    ax2.set_ylabel('ESI-MS signal intensity', color='b')

    ax.yaxis.label.set_color('g')
    ax.tick_params(axis='y', colors='g')

    ax2.yaxis.label.set_color('b')
    ax2.tick_params(axis='y', colors='b')

    
    if plt_time_limits and (icpms_max == None):
        ax.set_xlim(plt_time_limits[0], plt_time_limits[1])
        icpsub_lim = subset_icpdata(icp_data = icpdata, offset=0,timerange=[plt_time_limits[0],plt_time_limits[1]],element=heteroAtom)
        ymax =1.1*(max(icpsub_lim[heteroAtom]))
        ax.set_ylim(-0.005*ymax, ymax)
        ax2.set_ylim(-0.005*max(eics),max(eics)*1.1)

    elif plt_time_limits and (icpms_max != None):
        ax.set_xlim(plt_time_limits[0], plt_time_limits[1])
        ax.set_ylim(-0.005*icpms_max, icpms_max)
        ax2.set_ylim(-0.005*max(eics),max(eics)*1.1)
    
    elif plt_time_limits == None:
        ax.set_xlim([timestart, timestop])

    ax.text(0.05, 0.90, '$R^2$: %.3f' % cor_mz,transform=ax.transAxes, fontsize = 'small')

    if assigned and (pd.isna(mf) is False):    
        ax.text(0.05, 0.95, 'm/z: %.6f' %mz,transform=ax.transAxes, fontsize = 'small')
        ax.text(0.05, 0.85, 'Formula: %s' %(mf),transform=ax.transAxes, fontsize = 'small')
        ax.text(0.05,0.8, 'm/z error (ppm): %.3f' %mz_er,transform=ax.transAxes, fontsize = 'small')
        ax.text(0.05, 0.75, 'Confidence: %.3f' %cs,transform=ax.transAxes, fontsize = 'small')

    elif assigned and (pd.isna(mf) is True):    
        ax.text(0.05, 0.95, 'm/z: %.6f' %mz,transform=ax.transAxes, fontsize = 'small')
        ax.text(0.05, 0.85, 'Formula not assigned',transform=ax.transAxes, fontsize = 'small')

    ax.legend(handles=[line1, line2], loc = legendloc, frameon = False, fontsize = 'small')   

    ms_i = scan_range_data_i[0]
    print(scanrange)

    print(ms_i.mz_exp)

    ax_ms = fig.add_subplot(gs[ 1])  
    plot_ms(masspectrum = ms_i, srange=[scanrange, [timestart,timestop]], target_mz=mz, assignment=assignment_list, ax_ms=ax_ms)

    

    if show:
        #plt.tight_layout()
        plt.show()
    else:
        plt.subplots_adjust(bottom=bplt1, left = lplt1, right=rplt1, top=tplt1, hspace=hspace1, wspace =wspace1)
        dir = '/Users/christiandewey/Downloads/'
        plt.savefig(dir + title  + '.pdf',bbox_inches='tight')
    

def get_parMS_slist(scanlist_obj):

    scanlist = scanlist_obj.scanlist

    srange = scanlist_obj.srange
    scan_index = scanlist_obj.scan_index
    scans_per_chroma = scanlist_obj.scans_per_chroma
    scanlist = scanlist_obj.scanlist

    slist = []

    range_labels = []
    rmax = int(scan_index[scan_index[:,1] == scanlist[-1]][0,0] - scan_index[scan_index[:,1] == scanlist[0]][0,0])

    for r in srange: 
        
        if r == rmax:
            
            break

        if (r + scans_per_chroma) > rmax:
            
            final_r = rmax
        else:
            
            final_r = r+scans_per_chroma
        
        range_label = '%s-%s' % (scanlist[r], scanlist[final_r-1])

        range_labels.append(range_label)

        scanrange = scanlist[r:final_r]

        slist.append(scanrange)

    return slist, range_labels

def get_numscans_in_slist(slist):
    numscans = 0
    for i in slist:
        for j in i:
            numscans = numscans + 1
    return numscans

    

importlib.reload(parMS)


###
### wastewater sample
global ioncharge
ioncharge = 1

global bplt1
bplt1 = 0.01
global lplt1
lplt1 = 0.01
global rplt1
rplt1 = 0.99
global tplt1
tplt1 = 0.95
global wspace1
wspace1 = 0.5
global hspace1
hspace1 =0.3
global legendloc
legendloc = 'center left'
global mzrange 
mzrange = 4

import time
heteroAtom = '59Co'

MSParameters.mass_spectrum.threshold_method = 'signal_noise'
MSParameters.mass_spectrum.s2n_threshold = 6
MSParameters.ms_peak.peak_min_prominence_percent = 0.1 #0.1

timestart = 14.0
timestop = 15

offset =34.3 #-38 #Kansas soil, -27

cwindow = 20
elementDict_ww = {'C': (1,50), 'H':(4,100), 'O':(0,10), 'N':(0,5),'Co':(0,1)} # 'S':(0,4),


esifile_ww = '/Users/christiandewey/CoreMS/tests/tests_data/cobalt-wastewater/rmb_20220627_fp_CWD_wastewater_43.raw'
icpmsfile_ww = '/Users/christiandewey/CoreMS/tests/tests_data/cobalt-wastewater/cwd_220627_42_biozen_wastewater_50uL.csv'

icpdata_ww = pd.read_csv(icpmsfile_ww)
icpdata_ww.dropna(inplace=True)

esi_parser_ww = rawFileReader.ImportMassSpectraThermoMSFileReader(esifile_ww)

scanlist_obj_ww = ScanList( esi_parser_ww, icpdata_ww, chromaWindow=cwindow, timerange = [timestart,timestop] )


slist, range_labels = get_parMS_slist(scanlist_obj_ww)



#eics = extractEICs(scanlist_obj_ww)


eics_par = {}

for sl, rl, itl in zip(slist, range_labels, range(len(slist))):

    print('...extracting from scan range %s (%s / %s)' %(rl, itl+1, len(slist)))
    mass_spectrum = esi_parser_ww.get_average_mass_spectrum_by_scanlist(sl)

    eics_sl = parMS.run(esifile=esifile_ww, mz_exp=mass_spectrum.mz_exp) #nprocessors,

    eics_par[rl] = eics_sl


#scan_range_data_ww, range_labels_ww = process(elementDict = elementDict_ww, scanlist_obj=scanlist_obj_ww, eics=eics, assign_flag=True,returnLabels =True) #, specific_valence = ['Co',3]
scan_range_data_mp, range_labels_mp = process(elementDict = elementDict_ww, scanlist_obj=scanlist_obj_ww, eics=eics_par, assign_flag=True,returnLabels =True) #, specific_valence = ['Co',3]


threshold = 0.6

#best_ww, full_ww = getResults(esi_parser_ww, icpdata_ww, scan_range_data_ww, title = 'wastewater-ser',show=False, assigned = True, filter_corr = True, filter_hetero=True, plt_time_limits=[12,20], icpms_max=18000)
best_mp, full_mp = getResults(esi_parser_ww, icpdata_ww, scan_range_data_mp, title = 'wastewater-par',show=False, assigned = True, filter_corr = True, filter_hetero=True, plt_time_limits=[12,20], icpms_max=18000)












import gc

del heteroAtom
del esifile
del icpmsfile
del icpdata
del esi_parser
del scanlist_obj
del eics
del timestart
del timestop
del scan_range_data
del elementDict

gc.collect()


temp

mz_12c = temp['m/z'].values[3]
mz_13c = temp['m/z'].values[1]
mz_12c
rlabel = list(scan_range_data_ww.keys())[1]
scanrangeww = [2331,2411]

plot_single_mz(esi_parser_ww,icpdata_ww, scan_range_data_ww[rlabel],scanrangeww, mz=mz_13c, show = True, title = 'wastewater_470.241425',assigned = True, results = best_ww[1], plt_time_limits=[0,20], icpms_max=18000)

best_ww[1].to_csv('/Users/christiandewey/Downloads/bestww.csv')






#esifile = '/Users/christiandewey/Downloads/rmb_20220627_fp_1uM_stdmix_45.raw'



#esifile = '/Users/christiandewey/CoreMS/tests/tests_data/ftms/rmb_161221_kansas_h2o_2.raw'
#icpmsfile = '/Users/christiandewey/CoreMS/tests/tests_data/icpms/161220_soils_hypercarb_3_kansas_qH2O.csv'



heteroAtom = '127I'

MSParameters.mass_spectrum.threshold_method = 'signal_noise'
MSParameters.mass_spectrum.s2n_threshold = 3
MSParameters.ms_peak.peak_min_prominence_percent = 0.05 #0.1

timestart =44.6
timestop = 45.6

offset =-38 #

cwindow = 70





###
### 600 m depth
esifile = '/Users/christiandewey/CoreMS/tests/tests_data/marine-iodine/220822_CTD27_600m2.raw'
icpmsfile = '/Users/christiandewey/CoreMS/tests/tests_data/marine-iodine/CTD27_600m.csv'

icpdata = pd.read_csv(icpmsfile)
icpdata.dropna(inplace=True)

esi_parser = rawFileReader.ImportMassSpectraThermoMSFileReader(esifile)

scanlist_obj = ScanList( esi_parser, icpdata, chromaWindow=cwindow, timerange = [timestart,timestop] )
eics = extractEICs(scanlist_obj)

elementDict = {'C': (1,50), 'H':(4,100), 'O':(0,10), 'N':(0,10),'S':(0,4),'I':(0,1)}#'P':(0,4),

for r in scanlist_obj.srange:
    print(r)

scan_range_data, range_labels = process(elementDict = elementDict, scanlist_obj=scanlist_obj, eics=eics, assign_flag=True,returnLabels =True)

best600, full600 = getResults(esi_parser, icpdata, scan_range_data, title = '220822_CTD27_600m2',show=False, assigned = True, filter_corr = True, filter_hetero=False, plt_time_limits=[0,60])

ms = esi_parser.get_average_mass_spectrum_by_scanlist([8956,9146])
eic=esi_parser.get_eics(target_mzs=[mz],tic_data={},peak_detection=False,smooth=False)

###
### 1000 m depth
esifile = '/Users/christiandewey/CoreMS/tests/tests_data/marine-iodine/220822_CTD27_1000m2.raw'
icpmsfile = '/Users/christiandewey/CoreMS/tests/tests_data/marine-iodine/CTD27_1000m.csv'

icpdata1000 = pd.read_csv(icpmsfile)
icpdata1000.dropna(inplace=True)

esi_parser1000 = rawFileReader.ImportMassSpectraThermoMSFileReader(esifile)

scanlist_obj1000 = ScanList( esi_parser1000, icpdata1000, chromaWindow=cwindow, timerange = [timestart,timestop] )
eics1000 = extractEICs(scanlist_obj1000)

elementDict = {'C': (1,50), 'H':(4,100), 'O':(0,10), 'N':(0,10),'I':(0,1)}#'P':(0,4),
scan_range_data1000, range_labels1000 = process(elementDict = elementDict, scanlist_obj=scanlist_obj1000, eics=eics1000, assign_flag=True,returnLabels =True)

threshold = 0.4
best1000, full1000 = getResults(esi_parser1000, icpdata1000, scan_range_data1000,  title = '220822_CTD27_1000m2',show=False, assigned = True, filter_corr = True, filter_hetero=False, plt_time_limits=[0,60])





plot_single_mz(esi_parser, icpdata, mz=470.138058, show = False, title = '220822_CTD27_600m2_470.138058a',assigned = True, results = best600[0], plt_time_limits=[0,60], icpms_max=25000)


plot_single_mz(esi_parser1000, icpdata1000, mz=470.138423, show = False, title = '220822_CTD27_1000m2_470.138423b',assigned = True, results = best1000[0], plt_time_limits=[0,60], icpms_max=25000)

#plot_single_mz(esi_parser1000, mz=360.149755, show = True, title = '220822_CTD27_1000m2_360.236691',assigned = True, results = full1000, plt_time_limits=[0,60])





'''

import multiprocessing as mp
print("Number of processors: ", mp.cpu_count())

'''






len(scanlist_obj.scanlist)






