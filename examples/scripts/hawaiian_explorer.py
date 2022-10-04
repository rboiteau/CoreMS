
####Christian 09_30_2022 

import os
from matplotlib import style
import pandas as pd
import numpy as np

#Generate clusters of ms features across depth.
from sklearn.cluster import AgglomerativeClustering
import seaborn as sns
from matplotlib.colors import LogNorm, Normalize



import warnings

from examples.scripts.LCMS_assignAll import assign
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

import matplotlib.backends.backend_pdf

refmasslist = Path.cwd() / "tests/tests_data/ftms/nom_pos.ref"

#Set peak detection threshold method
MSParameters.mass_spectrum.threshold_method = 'signal_noise'
MSParameters.mass_spectrum.s2n_threshold = 3
MSParameters.ms_peak.peak_min_prominence_percent = 0.01

MSParameters.molecular_search.error_method = 'None'
MSParameters.molecular_search.min_ppm_error = -5
MSParameters.molecular_search.max_ppm_error = 5

MSParameters.molecular_search.url_database = None
MSParameters.molecular_search.min_dbe = -1
MSParameters.molecular_search.max_dbe = 20

MSParameters.molecular_search.usedAtoms['C'] = (1,50)
MSParameters.molecular_search.usedAtoms['H'] = (4,100)
MSParameters.molecular_search.usedAtoms['O'] = (1,20)
MSParameters.molecular_search.usedAtoms['N'] = (0,4)
#MSParameters.molecular_search.usedAtoms['S'] = (0,0)

MSParameters.molecular_search.isProtonated = True
MSParameters.molecular_search.isRadical = False
MSParameters.molecular_search.isAdduct = False

MSParameters.molecular_search.score_method = "prob_score"
MSParameters.molecular_search.output_score_method = "prob_score"

def get_b12_eic(parser):

    MSParameters.mass_spectrum.threshold_method = 'signal_noise'
    MSParameters.mass_spectrum.s2n_threshold =6
    MSParameters.ms_peak.peak_min_prominence_percent = 0.1

    tic = parser.get_tic(ms_type = 'MS')[0]

    mass_spectrum = parser.get_average_mass_spectrum_by_scanlist(tic.scans)

    mz_exp_np = np.array(mass_spectrum.mz_exp)

    sub = list(mz_exp_np[np.where((mz_exp_np > 678) & (mz_exp_np < 680))])

    b12_mass = sub
            
    '''for mi, mz in zip(range(len(sub)), sub):
        
        for i in range(mi+1,len(sub)):
            
            mz_i = sub[i]
            mz_diff = abs(mz_i - mz)

            if ( abs((mz_diff / 0.5014) - 1) * 10000) < 2 :
                
                b12_mass = mz 

                break
    '''
        
    b12EIC=parser.get_eics(target_mzs=b12_mass,tic_data={},peak_detection=False,smooth=False)

    return b12_mass, b12EIC
    
def qc_b12(parser,std_timerange, stdmass, linefmt=None):
    #QC control 

    
    ax=plt.subplot(111)
    ax.set(xlabel='Time (min)',ylabel='Intensity',title='Std EIC= %.4f' %(stdmass))

    EIC=parser.get_eics(target_mzs=[stdmass],tic_data={},peak_detection=False,smooth=False)
    df=pd.DataFrame({'EIC':EIC[0][stdmass].eic,'time':EIC[0][stdmass].time})
    df_sub=df[df['time'].between(std_timerange[0],std_timerange[1])]

    if linefmt:
        ax.plot(df_sub['time'],df_sub['EIC'],linefmt)
    else:
        ax.plot(df_sub['time'],df_sub['EIC'])
    
    area=sum(df_sub['EIC'])
    rt=df_sub.time[df_sub.EIC==df_sub.EIC.max()].max()

    return ax, area, rt

def run_b12_qc(MSfiles, stdmass, filelist, peakrange):

    plt.figure(figsize=(8,5))

    for file in MSfiles.keys():

        parser = MSfiles[file]['parser']

        if 'qH2O' in file:
            continue
        elif ('StdMix' in file) | ('stdmix' in file):
            ax, area, rt = qc_b12(parser,peakrange,stdmass,linefmt=':k')
        else:
            ax, area, rt = qc_b12(parser,peakrange,stdmass)

        MSfiles[file]['fig'] = ax
        MSfiles[file]['b12_area'] = area
        MSfiles[file]['b12_rt'] = rt
    
    ax.legend(labels = list(MSfiles.keys()), loc='center left', prop={'size': 6}, bbox_to_anchor=(1, 0.5))

    fig = plt.gcf()
    fig.tight_layout()
    fig.savefig('b12_eics.pdf',bbox_inches = 'tight')

    plt.show()
    
    areas = [MSfiles[f]['b12_area'] for f in filelist if ( ('qH2O' not in f) and ('StdMix' not in f) )]

    rts = [MSfiles[f]['b12_rt'] for f in filelist if ('qH2O' not in f) ]

    area_stdev = np.std(areas)
    area_mean = np.mean(areas)

    rt_stdev = np.std(rts)
    rt_mean = np.mean(rts)
    
    
    fig, (ax, axrt) = plt.subplots(1, 2, figsize=(6,4))

    ax.boxplot(areas, labels=['positive mode'])
    ax.set_ylabel('b12 peak area')
    ax.text(x = .05  , y = 0.95 , s = '%.1f%% std dev' %(area_stdev/area_mean*100),transform=ax.transAxes, fontsize = 'x-small')
    axrt.boxplot(rts, labels=['positive mode'])
    axrt.set_ylabel('b12 r.t. (min)')
    axrt.text(x = .05  , y = 0.95 , s = '%.1f%% std dev' %(rt_stdev/rt_mean*100),transform=axrt.transAxes, fontsize = 'x-small')
    fig.tight_layout()
    fig.savefig('b12_qcs.pdf')
    #plt.tight_layout()
    plt.show()

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

def plot_reqd_resolving_power(masterresults):

    #Make plots showing required resolving power.
    only_samples = {k: masterresults[k] for k in masterresults.keys() if ( 'stdmix' not in k) and ('blank' not in k)  }
    masterresults_df=pd.concat(only_samples.values())

    diff_summary=[]
    res_summary=[]
    for time in masterresults_df['Time'].unique():
        result=masterresults_df[masterresults_df['Time']==time]

        for file in masterresults_df['file'].unique():
            result_sub=result[result['file']==file]
            result_sub=result_sub['m/z'].sort_values(ascending=True)
            #mzdiff=result_sub['m/z'].sort_values(ascending=True).diff().iloc[1:]
            #mzdiff_res=result_sub['m/z'].iloc[1:]/mzdiff
            differences=result_sub.diff()


            #Resolve from peaks on either side:
            mzdiff=pd.DataFrame({'left':differences[1:-1].to_list(),'right':differences[2:].to_list(),'mz':result_sub[1:-1].to_list()})
            mzdiff_min=abs(mzdiff.min(axis=1))
            mzdiff_res=mzdiff['mz']/(mzdiff_min/2)

            diff_summary.extend(mzdiff_min.tolist())
            res_summary.extend(mzdiff_res[mzdiff_res<1e6].tolist())
            #res_summary.extend(mzdiff_res.tolist())

    #print(diff_summary)

    res=list(range(10,int(max(res_summary)),1000))

    count=[]
    for i in res:
        count.append(len([element for element in res_summary if element<i])/len(res_summary))

    #print(len(diff_summary))
    #print(len([element for element in diff_summary if element<0.01]))
    fig, ax = plt.subplots()
    ax.plot(res,count)
    for i in [0.8,0.9,0.95,0.99]:
        current=abs(np.array(count)-i).tolist()
        print(res[current.index(min(current))])
        ax.axvline(res[current.index(min(current))],color='black')
    ax.set_xlabel('Resolving power (m/z / dm/z)')
    ax.set_ylabel('Fraction of peaks resolved')

    #plt.hist(diff_summary, bins=50,range=[0,.03])
    #plt.xlabel('m/z diff (Da)')
    #plt.ylabel('frequency')
    return ax

def plot_results(masterresults):

    only_samples = {k: masterresults[k] for k in masterresults.keys() if ( 'stdmix' not in k) }
    masterresults_df=pd.concat(only_samples.values())

    allresults=masterresults_df[(masterresults_df['m/z']<800) & (masterresults_df['S/N']>3)]

    allresults['N']=allresults['N'].fillna(0)
    allresults['O']=allresults['O'].fillna(0)
    #allresults['S']=allresults['S'].fillna(0)
    #allresults['P']=allresults['P'].fillna(0)
    #allresults['Na']=allresults['Na'].fillna(0)

    allresults['mol_class']='Unassigned'
    allresults['mol_class'][allresults['C']>0]='CHO'
    allresults['mol_class'][(allresults['C']>0) & (allresults['N']>0.5)]='CHON'
    #allresults['mol_class'][(allresults['C']>0) & (allresults['S']>0.5)]='CHOS'
    #allresults['mol_class'][(allresults['C']>0) & (allresults['P']>0.5)]='CHOP'
    #allresults['mol_class'][(allresults['C']>0) & (allresults['Na']>0.5)]='CHONa'
    #allresults['mol_class'][(allresults['C']>0) & (allresults['S']>0.5) & (allresults['N']>0.5)]='CHONS'
    #allresults['mol_class'][(allresults['C']>0) & (allresults['P']>0.5) & (allresults['N']>0.5)]='CHONP'
    #allresults['mol_class'][(allresults['C']>0) & (allresults['Na']>0.5) & (allresults['N']>0.5)]='CHONNa'
    #allresults['mol_class'][(allresults['C']>0) & (allresults['P']>0.5) & (allresults['Na']>0.5) & (allresults['N']>0.5)]='CHONPNa'


    results= allresults[allresults['Is Isotopologue']==0]
    results['O/C']=results['O']/results['C']
    results['H/C']=results['H']/results['C']
    results['N/C']=results['N']/results['C']


    colors = {'CHO':'red', 'CHON':'green', 'CHOS':'blue', 'CHONS':'yellow', 'CHONP':'black', 'CHONNa':'cyan','CHONPNa':'pink','CHONa':'aqua','CHOP':'gray'}

    grouped=results.groupby('mol_class')
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1,figsize=(8,12))
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
    for time in allresults['Time'].unique():
        current={}
        current['Time']=time
        for mol_class in sorted(allresults['mol_class'].unique()):
            current[mol_class]=len(allresults[(allresults['mol_class']==mol_class) & (allresults['Time']==time)])
        assign_summary.append(current)
    df=pd.DataFrame(assign_summary)

    ax4 = df.plot.bar(x='Time',y=df.columns[1:],stacked=True,ylabel='Peaks',ax=ax4)
    ax4.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    return [ax1, ax2, ax3, ax4], allresults

def plot_unique_features(allresults):

    #Create a list of all unique features and describe their intensity. 
    print(len(allresults))
    #define a list of unique features (time, formula) with 'areas' determined for each sample. There may be a slight bug that causes the unique list to grow...
    uniquelist=[]
    for time in allresults.Time.unique():
        current=allresults[allresults.Time==time]
        current=current.sort_values(by=['m/z Error (ppm)'],ascending=True)
        currentunique=current.drop_duplicates(subset=['Molecular Formula'])
        currentunique=currentunique[currentunique['C']>1]
        currentunique=currentunique.set_index(['Molecular Formula'],drop=False)
        for file in allresults['file'].unique():
            current_file=current[current['file']==file].drop_duplicates(subset=['Molecular Formula'])
            current_file=current_file.rename(columns={'Peak Height':file})
            current_file=current_file.set_index(['Molecular Formula'],drop=False)
            currentunique=currentunique.join(current_file[file])
        uniquelist.append(currentunique)

    uniqueresults=pd.concat(uniquelist,ignore_index=True)
    uniqueresults['N/C']=uniqueresults['N']/uniqueresults['C']
    uniqueresults['blank']=uniqueresults['rmb_20220627_fp_blank1_22.raw']/uniqueresults['Peak Height']
    #uniqueresults=uniqueresults[uniqueresults['blank']<0.95]


    print(len(uniqueresults))


    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1,figsize = (8,12))
    sns.violinplot(x="Time", y="O/C", hue='mol_class', data=uniqueresults, ax=ax1)
    sns.violinplot(x="Time", y="N/C", hue='mol_class', data=uniqueresults, ax=ax2)
    ax.set(xlabel='Time (min)')
    sns.scatterplot(x='m/z',y='m/z Error (ppm)',hue='N',data=uniqueresults,ax=ax3,s=10)
    sns.scatterplot(x='m/z',y='Resolving Power',hue='mol_class',data=uniqueresults,ax=ax4)
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax4.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    return [ax1,ax2,ax3,ax4], uniqueresults

def run_cluster_analysis(uniqueresults, masterresults, sort_order):

    #Generate clusters of ms features across depth.
    from sklearn.cluster import AgglomerativeClustering
    import seaborn as sns
    from matplotlib.colors import LogNorm, Normalize


    #Blank filtering:

    only_samples = {k: masterresults[k] for k in masterresults.keys() if   ('stdmix' not in k) and ('blank' not in k)}
    masterresults_df=pd.concat(only_samples.values())

    #clustermethod='average'
    clustermethod='ward'

    abundances=uniqueresults[masterresults_df['file'].unique()].fillna(0)

    df=abundances.mean(axis=1)
    df_std=abundances.std(axis=1)

    p_list=[]
    for ind in abundances.index:
        p=max(abs(abundances.loc[ind]-df[ind]))/df_std[ind]
        if len(abundances.loc[ind][abundances.loc[ind]>1])<2:
            p=0
        p_list.append(p)

    #abundances['p']=p_list
    #abundances=abundances[abundances['p']>0.1]
    #abundances=abundances.sort_values(by='p',ascending=False)

    uniqueresults['p']=p_list

    #plt.hist(p_list,bins=100,range=[min(p_list),max(p_list)])

    results_clustered=uniqueresults[uniqueresults['p']>0.1]
    norm_abundances=results_clustered[masterresults_df['file'].unique()].fillna(0)
    norm_abundances=norm_abundances.div(norm_abundances.max(axis=1),axis=0)


    cluster = AgglomerativeClustering(n_clusters=3,affinity='euclidean',linkage=clustermethod)
    cluster.fit_predict(norm_abundances)

    results_clustered['cluster']=cluster.labels_

    #results_clustered.fillna(0).to_csv(file_location+'clustered_results.csv')



    clusterplot=norm_abundances
    #clusterplot=clusterplot.drop(['p'],axis=1)
    clusterplot=clusterplot.transpose()
    clusterplot['sort_order']=sort_order
    clusterplot=clusterplot.sort_values(by='sort_order')
    clusterplot=clusterplot.drop(['sort_order'],axis=1)
    p1=sns.clustermap(clusterplot,row_cluster=False,cmap='mako',method=clustermethod)

    plt.show()







###### Set file folder and THERMO RAW file name here:
file_location='/Users/christiandewey/Downloads/Hawaiian_soils/Thermo_RAW/pos/'
#file_location='/Users/christiandewey/Library/CloudStorage/Box-Box/Boiteau Lab/Mass Spec Data/21T at NHMFL/2021_August_OC2102A/OC2012A/'
filelist=os.listdir(file_location)
os.chdir(file_location)

MSfiles={}
for file in filelist:
    if '.raw' in file:
        parser = rawFileReader.ImportMassSpectraThermoMSFileReader(file_location+'/'+file)
        MSfiles[file]={'parser': parser}



##### get EIC of b12 std 



std_parser = MSfiles['rmb_20220627_fp_1uM_stdmix_45.raw']['parser']

stdmass, stdEIC = get_b12_eic(std_parser)

fig,ax = plt.subplots()

max_std_intensity = 0
dominant_mz = None
for mz in stdmass:
    ax.plot(stdEIC[0][mz].time, stdEIC[0][mz].eic, label='%.4f ' %mz)
    ax.legend(loc='best', frameon=False)

    if max_std_intensity < max(stdEIC[0][mz].eic):

        max_std_intensity = max(stdEIC[0][mz].eic)

        dominant_mz = mz

plt.show()

max_std_intensity
dominant_mz



##### run b12 qc 
raw_files = [f for f in filelist if '.raw' in f]

run_b12_qc(MSfiles,dominant_mz,raw_files,[8,8.5])


#### assign formula 

timerange = [2,25]
interval = 1


master_results = {}
for file in MSfiles:
    results = assign_formula(MSfiles[file]['parser'], interval, timerange, refmasslist=refmasslist)
    results['file'] = file 
    MSfiles[file]['results'] = results
    master_results[file] = results

ax = plot_reqd_resolving_power(master_results)

plt.tight_layout()
plt.savefig('%sresolving_power.pdf' %(file_location))



axs, allresults = plot_results(master_results)

plt.tight_layout()
#plt.show()
plt.savefig('%sall_results.pdf' %(file_location))



u_axs, unique_results = plot_unique_features(allresults)

plt.tight_layout()
plt.show()



sort_order=[]
for file in filelist:
    if ('fp_p' in file): # or ('blank' in file):
        d=file.split('_')[3]
        sort_order.append(d)

run_cluster_analysis(unique_results,master_results,sort_order)