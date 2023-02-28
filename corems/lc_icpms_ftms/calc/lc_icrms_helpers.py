
####Christian 12_13_2022 

import os
from matplotlib import style
import pandas as pd
import numpy as np
import tqdm

#Generate clusters of ms features across depth.
from sklearn.cluster import AgglomerativeClustering
import seaborn as sns
from matplotlib.colors import LogNorm, Normalize

import warnings

warnings.filterwarnings("ignore")
import sys
sys.path.append("./")
from pathlib import Path
import matplotlib.pyplot as plt

from corems.mass_spectra.input import rawFileReader
from corems.molecular_id.factory.classification import HeteroatomsClassification, Labels
from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.encapsulation.factory.parameters import MSParameters
from corems.mass_spectrum.calc.Calibration import MzDomainCalibration


"""
Christian Dewey 
December, 2022

Helper functions for processing 21T data for m/z windowing project

Data collected from NHMFL in November, 2022
"""

def getParser(file):
    parser = rawFileReader.ImportMassSpectraThermoMSFileReader(file)
    return parser


   


def plot_ms(df1, start_mz, end_mz,mfexclude = [],tstart=None, df2=None,df3=None, assignment= None, ax_ms=None, lbls=None, norm=False, labs=False, colors=None):   
    
    if ax_ms == None:
        _, ax = plt.subplots()
    
    else:
        ax = ax_ms

    if colors is not None:
        cols = colors
    else:
        cols = ['C0', 'C1', 'C2']
    mzrange= end_mz - start_mz

    if tstart != None:
        ms_t_int=df1[df1['Time'] == tstart]
    else:
        ms_t_int=df1

    ms_df = ms_t_int[((ms_t_int['Calibrated m/z']-start_mz)<mzrange) & ((ms_t_int['Calibrated m/z']-start_mz)>0)]


    if norm:

        pltcol = 'Normalized Peak Height'

    else:
        pltcol = 'S/N'

    
    if lbls is not None:
        labels = lbls
    else:
        labels = [None, None, None]

 
    _, stemlines1, _ =ax.stem('Calibrated m/z',pltcol,data=ms_df[~ms_df['Molecular Formula'].isin(mfexclude)],  markerfmt=' ', basefmt=' ', linefmt=cols[0], label = labels[0])
    
    if (df2 is not None) and (len(df2['Molecular Formula'])>0):
        pltdf2 = True
        if tstart != None:
            ms_t_int2=df2[df2['Time'] == tstart]
        else:
            ms_t_int2=df2

        ms_df2 = ms_t_int2[(abs(ms_t_int2['Calibrated m/z']-start_mz)<mzrange)& ((ms_t_int2['Calibrated m/z']-start_mz)>0)]
        
        _, stemlines2, _ =ax.stem('Calibrated m/z',pltcol,data=ms_df2[~ms_df2['Molecular Formula'].isin(mfexclude)],  markerfmt=' ', basefmt=' ', linefmt=cols[1],label = labels[1])

         #markerline, stemlines, baseline = plt.stem(x, y)
    else:
        pltdf2 = False 

    if (df3 is not None) and (len(df3['Molecular Formula'])>0):
        pltdf3 = True
        if tstart != None:
            ms_t_int3=df3[df3['Time'] == tstart]
        else:
            ms_t_int3=df3

        ms_df3 = ms_t_int3[(abs(ms_t_int3['Calibrated m/z']-start_mz)<mzrange)& ((ms_t_int3['Calibrated m/z']-start_mz)>0)]

        
        _, stemlines3, _ =ax.stem('Calibrated m/z',pltcol,data=ms_df3[~ms_df3['Molecular Formula'].isin(mfexclude)],  markerfmt=' ', basefmt=' ', linefmt=cols[2], label = labels[2])
    else:
        pltdf3 = False
    #if pltdf3 is True:
    #    ax.set_ylim(0, max([maxdf1, maxdf2, maxdf3]) * 1.1)
    #elif pltdf2 is True:
    #    ax.set_ylim(0, max([maxdf1, maxdf2]) * 1.1)
    #else: 
    #    ax.set_ylim(0, maxdf1 * 1.1)

    ax.set_xlim(left = start_mz - mzrange*0.1, right = start_mz + mzrange + mzrange*0.1) 

    if labs:
        for mzr,peakr,mf,er in zip(ms_df['Calibrated m/z'], ms_df[pltcol], ms_df['Molecular Formula'],  ms_df['m/z Error (ppm)']):

            #if (mzr- target_mz)  == 0:
            #    mz_text = ' m/z\n%.4f' % (mzr)
            #else:
            #    mz_text = r'$\Delta$' + ' m/z\n%.4f' % (mzr- target_mz)

            mz_text = 'm/z %.4f\n%s\n%.3f ppm' % (mzr,mf,er)
            ax.text(mzr, peakr + 0.02 *max(ms_df[pltcol]), mz_text, ha = 'center', fontsize = 'xx-small', weight = 'bold', color=cols[0])

        if df2 is not None:

            for mzr,peakr,mf, er in zip(ms_df2['Calibrated m/z'], ms_df2[pltcol], ms_df2['Molecular Formula'], ms_df2['m/z Error (ppm)']):

            #if (mzr- target_mz)  == 0:
            #    mz_text = ' m/z\n%.4f' % (mzr)
            #else:
            #    mz_text = r'$\Delta$' + ' m/z\n%.4f' % (mzr- target_mz)

                mz_text = 'm/z %.4f\n%s\n%.3f ppm' % (mzr,mf,er)
                ax.text(mzr, peakr + 0.02 *max(ms_df2[pltcol]), mz_text, ha = 'center', fontsize = 'xx-small', weight = 'bold', color = cols[1])

        if df3 is not None:

            for mzr,peakr,mf, er in zip(ms_df3['Calibrated m/z'], ms_df3[pltcol], ms_df3['Molecular Formula'], ms_df3['m/z Error (ppm)']):

            #if (mzr- target_mz)  == 0:
            #    mz_text = ' m/z\n%.4f' % (mzr)
            #else:
            #    mz_text = r'$\Delta$' + ' m/z\n%.4f' % (mzr- target_mz)

                mz_text = 'm/z %.4f\n%s\n%.3f ppm' % (mzr,mf,er)
                ax.text(mzr, peakr + 0.02 *max(ms_df3[pltcol]), mz_text, ha = 'center', fontsize = 'xx-small', weight = 'bold', color = cols[2])

   # theor_mz=pattern.mdiff+result['mass']
   # theor_int=pattern.ratio*result['abundance']
   # ax.stem(theor_mz,theor_int, basefmt=' ',linefmt='gray')

   # for isotope in pattern.isotope[pattern.requirement=='Y']:
   #     ax.stem('mz','intense',data=result[isotope],  markerfmt=' ', basefmt=' ',linefmt='red')
    ax.set_xlim(start_mz, end_mz)
    
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
        ax.set(xlabel='Calibrated m/z',ylabel='Normalized Intensity')
    else: 
        ax.set(xlabel='Calibrated m/z',ylabel='S/N')
    #ax.set_title('%.2f' %timerange[0] + ' to %.2f' %timerange[1] +' min', fontsize = 'medium')
    ax.legend(bbox_to_anchor = (1.00, 0.5), frameon =False, loc = 'center left')
    ax.axhline(y=0.0, color='black')
    plt.setp(stemlines1,'color', cols[0], 'linewidth', 2)
    if pltdf2 is True:
        plt.setp(stemlines2, 'color', cols[1],'linewidth', 1.75, 'linestyle', '--')
    if pltdf3 is True:
        plt.setp(stemlines3, 'color', cols[2],'linewidth', 2)
    plt.tight_layout()
    if ax_ms != None:
        return ax


def filterMzRange(results, mz_range):

    mz_i = mz_range[0]
    mz_f = mz_range[1]

    sub = results[(results['m/z'] >= mz_i) & (results['m/z'] <= mz_f)]

    return sub


def pltMZerror(results, bins=50):
    ## assignment error distribution
    #_, ax = plt.subplots()
    #for mol_class in sorted(results['mol_class'].unique()):

        #counts, bins = np.histogram(np.asarray(results[results['mol_class']==mol_class]['m/z Error (ppm)']),bins = bins)

        #ax.plot(bins[:-1], counts, label = mol_class)

    sns.kdeplot(data=results, x="m/z Error (ppm)", hue="mol_class")

    #ax.set_xlim(-0.5,0.5)
    #ax.legend(frameon=False)
    #ax.set_xlabel('m/z assignment error (ppm)')
    #ax.set_ylabel('Density')

    #return ax


def pltMZerror_pts(results):
    ## assignment error distribution
    _, ax = plt.subplots()
    for mol_class in sorted(results['mol_class'].unique()):

        ax.scatter(results[results['mol_class']==mol_class]['Calibrated m/z'], results[results['mol_class']==mol_class]['m/z Error (ppm)'], label = mol_class)

    ax.set_xlim(200,1200)
    ax.legend(frameon=False,bbox_to_anchor=(1.0, 0.5))
    ax.set_xlabel('calibrated m/z')
    ax.set_ylabel('assignment error (ppm)')

    return ax


def assign_mol_class(complete_results_df,molclasses,sn_threshold=3,mz_threshold=800):

    all_results=complete_results_df[(complete_results_df['m/z']<mz_threshold) & (complete_results_df['S/N']>sn_threshold)]

    i = 0 # iterable 
    p = 0 # max len holder
    j = 0 # index of max len mol class

    for i in range(0,len(molclasses)):

        mc = molclasses[i]
        
        if (len(mc) > p) and (mc != 'Unassigned'):
            
            p = len(mc)
            
            j = i

    all_elements = get_elements(molclasses[j])

    all_results['ID'] = range(0,len(all_results))
    
    times = all_results['Time'].unique()


    holder = []
    sizenp = 0
    for t in times:
        print('\ntime:',t)
        time_average = all_results[all_results['Time'] == t]

        sizenp = sizenp + len(time_average)

        print('unassigned: ', len(time_average[time_average['Molecular Formula'].isna()]))
        print('assigned: ', len(time_average[~time_average['Molecular Formula'].isna()]))

        for m in molclasses:

            if m != 'Unassigned':
                elements = get_elements(m)
                
                sub = get_molclass_subset(elements, all_elements,time_average[~time_average['Molecular Formula'].isna()]) 
                sub['mol_class'] = m    

 
            elif m == 'Unassigned':
                sub = time_average[time_average['Molecular Formula'].isna()] 
                sub['mol_class'] = m

            print('\t%s: %s' %(m,len(sub)))

            holder.append(sub)


    results = pd.concat(holder)
    
    return results 


def _get_mol_classes(add, base):
    
    new = []
    remain = []
    
    new.append(base)

    for i in range(len(add)):
        
        new.append(base + add[i])

        new2 = []

        remain = add[i+1:]

        for j in remain:

            new2 = add[i] + j

            new.append(base + new2)

    return(new)


def get_mol_class(add):
    base = 'CHO'
    molclasses = []
    for i in range(len(add)):
        e = add[i]
        temp = _get_mol_classes(add[i:], base = base)
        base = base+e
        molclasses = molclasses + temp


    output = []
    for x in molclasses:
        if x not in output:
            output.append(x)

    output.append('Unassigned')
    return output


def get_elements(molclass):

    import re
    
    elements = [] 

    elements = re.findall('[A-Z][^A-Z]*', molclass)
    
    return elements

def get_molclass_subset(included_elements, all_elements, all_results):
    
    tdf = all_results

    for e in all_elements:
        tdf[e].fillna(0, inplace = True)

    excluded_elements = [e for e in all_elements if e not in included_elements]
    
    for e in included_elements:
        
        tdf = tdf[tdf[e]>0]
        for j in excluded_elements:
            tdf = tdf[tdf[j]==0]

    return tdf


def get_ratios(results):
    #results['H/C'] = 0
    ##results['O/C'] = 0
    #results['N/C'] = 0

    results[results['Is Isotopologue']==0]['O/C'] = results['O'] / results['C']
    results[results['Is Isotopologue']==0]['H/C'] = results['H'] / results['C']
    results[results['Is Isotopologue']==0]['N/C'] = results['N'] / results['C']

    return results 

def add_mzwindow_col(df):    

    df['m/z window'] = df.index
    df['Window Size (m/z)'] = df.index
    
    for file, r in zip(df['file'], range(len(df['file']))):

        if 'StdMix' not in file:

            if ('400_500' in file) or ('400-500' in file):

                df['m/z window'].iloc[r] = '400-500 m/z'
                df['Window Size (m/z)'].iloc[r] = "100"

            elif ('500_600' in file) or ('500-600' in file):

                df['m/z window'].iloc[r] = '500-600 m/z'
                df['Window Size (m/z)'].iloc[r] = "100"
        
            elif ('600_700' in file) or ('600-700' in file):

                df['m/z window'].iloc[r] = '600-700 m/z'
                df['Window Size (m/z)'].iloc[r] = "100"

            elif ('700_800' in file) or ('700-800' in file):

                df['m/z window'].iloc[r] = '700-800 m/z'
                df['Window Size (m/z)'].iloc[r] = "100"
            
            elif ('300_500' in file) or ('300-500' in file):

                df['m/z window'].iloc[r] = '300-500 m/z'
                df['Window Size (m/z)'].iloc[r] = "200"

            
            elif ('400_600' in file) or ('400-600' in file):

                df['m/z window'].iloc[r] = '400-600 m/z'
                df['Window Size (m/z)'].iloc[r] = "200"
            
            elif ('600_800' in file) or ('600-800' in file):

                df['m/z window'].iloc[r] = '600-800 m/z'
                df['Window Size (m/z)'].iloc[r] = "200"
                
            elif 'full' in file:

                df['m/z window'].iloc[r] = '200-1200 m/z'
                df['Window Size (m/z)'].iloc[r] = "1000"
        
        else:
            df['m/z window'].iloc[r] = '200-1200 m/z'
            df['Window Size (m/z)'].iloc[r] = "1000"

    return df 




def getUniqueFeatures(df):    
    '''
    Notes by Christian Dewey, 23-Feb-27

    OVERVIEW:
    
        - subset m/z assignments by time bin
        - sort subset by m/z error (ascending)
        - remove duplicate m.f. assignments, preserving the assignment with lowest m/z error (first hit on sorted subset), save as 'currentunique'
        - set the index of 'currentunique' to the 'Molecular Formula' column
        - for each file in the raw file list:
            - save subset the sorted m/z assignment list (with all files and duplicates) by file name as 'current_file'
            - remove duplicate formulae from file subset copies
            - rename 'Peak Height' to filename for each file
            - set 'Molecular Formula' to index for copy with renamed 'Peak Height' col
            - join the renamed 'Peak Height' col to 'currentunique'

    RETURNS:    Pandas dataframe 
                Contains unique m.f. assignments in each time bin, and intenisty of unique m.f. feature in each file within each time bin

    
    '''
    uniquelist=[]
    for time in df.Time.unique():
        current=df[df.Time==time]
        current=current.sort_values(by=['m/z Error (ppm)'],ascending=True)
        currentunique=current.drop_duplicates(subset=['Molecular Formula'])
        currentunique=currentunique[currentunique['C']>1]
        currentunique=currentunique.set_index(['Molecular Formula'],drop=False)
        for file in df['file'].unique():
            current_file=current[current['file']==file].drop_duplicates(subset=['Molecular Formula'])
            current_file=current_file.rename(columns={'Peak Height':file})
            current_file=current_file.set_index(['Molecular Formula'],drop=False)
            currentunique=currentunique.join(current_file[file])
        for mzw in df['Window Size (m/z)'].unique():
            current_file=current[current['Window Size (m/z)']==mzw].drop_duplicates(subset=['Molecular Formula'])
            wlbl = mzw + ' m/z window'
            current_file=current_file.rename(columns={'Peak Height':wlbl})
            current_file=current_file.set_index(['Molecular Formula'],drop=False)
            currentunique=currentunique.join(current_file[wlbl])
        uniquelist.append(currentunique)

    unique_results=pd.concat(uniquelist,ignore_index=True)
    unique_results['N/C']=unique_results['N']/unique_results['C']

    return unique_results



def plotUnique1(df,ps=50,includeBlanks=False, xlim = [400,600]):
    xmin = xlim[0]
    xmax = xlim[1]
    if includeBlanks != True:
        mask = ~df['file'].str.contains('qh2o', case=False, na=False)
        df=df[mask]
    else:
        df=df
    fig, ((ax2, ax4),(ax3, ax1)) = plt.subplots(2,2)
    sns.scatterplot(x='m/z',y='m/z Error (ppm)',data=df[df['Cu']>0],hue='m/z window', s=ps*4, ax=ax1)
    plt.legend(frameon=False)

    sns.scatterplot(x='m/z',y='m/z Error (ppm)', hue='m/z window', data=df, ax=ax2,s=ps)
    plt.legend(frameon=False)

    ax1.set_xlim(xmin,xmax)
    ax2.set_xlim(xmin,xmax)
    sns.scatterplot(x='m/z',y='m/z Error (ppm)',hue='N',data=df,ax=ax3,s=ps)
    plt.legend(frameon=False)

    sns.scatterplot(x='m/z',y='Resolving Power',hue='mol_class',data=df,ax=ax4)
    
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax1.set_title('Assignment error, Cu features')
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax2.set_title('Overall assignment error')
    ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax3.set_title('Assignment error, N features')
    ax4.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax4.set_title('Resolving power v. m/z')

    #fig.savefig('unique_results.pdf', bbox_to_inches='tight')

    return fig


def plotUnique(df,ps=50,includeBlanks=False, xlim = [400,600]):

    import matplotlib as mpl
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap
    cmap = sns.color_palette(["darkorange", "firebrick", "lightseagreen", "black","red",'yellow','blue','magenta','green','salmon'])
    #cmap = sns.color_palette("cubehelix_r")
    xmin = xlim[0]
    xmax = xlim[1]

    if includeBlanks != True:
        mask = ~df['file'].str.contains('qh2o', case=False, na=False)
        df=df[mask]
    else:
        df=df
    fig, (ax2, ax4, ax3) = plt.subplots(3,1,figsize = (12,12))
    sns.scatterplot(x='m/z',y='m/z Error (ppm)', hue='m/z window', data=df.sort_values(by='m/z window'), ax=ax2,s=ps)

    df['N'] = df['N'].fillna(0)
    df = df.astype({"N":'int'}) 

    print(np.sort(df['N'].unique()))
    sns.scatterplot(x='m/z',y='m/z Error (ppm)',hue='N',data=df[df['N']>0].sort_values(by='N'),palette = cmap,ax=ax3,s=ps)

    sns.scatterplot(x='m/z',y='Resolving Power',hue='mol_class',data=df,ax=ax4)

    
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon = False, title = 'm/z Window')
    ax2.set_title('Overall assignment error')



    ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon = False, title = 'N Atoms')
    ax3.set_title('Assignment error, N features')
    ax4.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon = False, title = 'Mol. Class')
    ax4.set_title('Resolving power v. m/z')

    #fig.savefig('unique_results.pdf', bbox_to_inches='tight')

    return fig



def addRepCol(data_df):

    data_df['Rep'] = data_df.index


    for file in data_df['file'].unique():

        print(file)

        if ('rep2' in file) or ('_02.' in file):

            temp = data_df[data_df['file'] == file]
            temp['Rep'] = 2
            data_df[data_df['file'] == file] = temp


        else:

            temp = data_df[data_df['file'] == file]
            temp['Rep'] = 1
            data_df[data_df['file'] == file] = temp

    print(data_df['Rep'].unique())


    return data_df 


def add_mz_window_colsl(data_df):

    data_df['Rep'] = data_df.index


    for file in data_df['file'].unique():

        print(file)

        if ('rep2' in file) or ('_02.' in file):

            temp = data_df[data_df['file'] == file]
            temp['Rep'] = 2
            data_df[data_df['file'] == file] = temp


        else:

            temp = data_df[data_df['file'] == file]
            temp['Rep'] = 1
            data_df[data_df['file'] == file] = temp

    print(data_df['Rep'].unique())


    return data_df 



def blankSubtract(df, blnkthresh = 0.8):
    # must be performed on df with unique assignments 
    holder = []
    for file in df['file'].unique():
        
        sub = df[df['file'] == file]

        blkf = sub['blank file'].iloc[0]

        sub[sub[file]== np.nan] = 0  # each file column contains intensities of feature; if feature is not present, nan assigned; need to convert to zero for blank subtract

        nom = sub[file]
        den = sub[blkf]

        nom = nom.replace(np.nan,0)  # features not present in sample
        den = den.replace(np.nan,1)  # features not present in blank

        if file != blkf:
            nom = nom
        elif file == blkf:
            nom = nom * (blnkthresh * 0.8)  # multiplication enables removal of these features from blank file 

        sub['blank subtract'] = nom/den  # ratio of intensities of features in sample and blank

        holder.append(sub)

    df_end = pd.concat(holder)

    df_end = df_end[df_end['blank subtract'] > blnkthresh]  # only select features that do not appear in blanks

    return df_end

def intersection(lst1, lst2):
    print(type(lst1))
    # Use of hybrid method
    temp = set(lst2)
    lst3 = [value for value in lst1 if value in temp]
    return lst3


def exclusion(ser1, ser2):
    # takes pd Series (e.g. DataFrame columns); returns list
    return list(set(list(ser1)) ^ set(list(ser2)))

def repCombine(df):

    for file in df['file'].unique():

        df[df[file] == np.nan] = 0

        if 'rep2' not in file:

            if '.raw' in file:

                rep2file = file.split('.')[0]+'_rep2.raw'

            else:

                rep2file = file + '_rep2'
            
            avfile = file + '_av'
            
            df[avfile] = (df[file] + df[rep2file]) / 2

    return df


def normMS(df,fulldf):

    max_i = max(fulldf['Peak Height'].values)

    df['Normalized Peak Height'] = df['Peak Height'] / max_i

    return df



def assign_formula(parser, interval, timerange, refmasslist=None,corder=2,charge=1, cal_ppm_threshold=(-1,1)):
    #Function to build formula assignment lists
    #Retrieve TIC for MS1 scans over the time range between 'timestart' and 'timestop' 
    tic=parser.get_tic(ms_type='MS')[0]
    tic_df=pd.DataFrame({'time': tic.time,'scan': tic.scans})

    times=list(range(timerange[0],timerange[1],interval))

    results=[]
    
    for timestart in times:
        print('timestart: %s' %timestart )
        scans=tic_df[tic_df.time.between(timestart,timestart+interval)].scan.tolist()

        mass_spectrum = parser.get_average_mass_spectrum_by_scanlist(scans)    
        mass_spectrum.molecular_search_settings.ion_charge = charge
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
                                                        calib_ppm_error_threshold=cal_ppm_threshold,
                                                        calib_snr_threshold=3)

            calfn.recalibrate_mass_spectrum(mass_spectrum, imzmeas, mzrefs, order=corder)


        SearchMolecularFormulas(mass_spectrum, first_hit=False).run_worker_mass_spectrum()

        mass_spectrum.percentile_assigned(report_error=True)

        assignments=mass_spectrum.to_dataframe()

        assignments['Time']=timestart

        results.append(assignments)
    
    results=pd.concat(results,ignore_index=True)

    return(results)    


 