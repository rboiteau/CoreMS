# Import the os module
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp
from scipy.stats import ranksums
import sys
sys.path.append("./")

from corems.mass_spectra.input import rawFileReader

### Set data directory and file names here

data_dir = '/CoreMS/usrdata/'
sample_list_name='BATS_sample_list.csv' #Sample list must contain column with header 'File'
featurelist_file='BATS_featurelist.csv'
clusteredlist_file='BATS_featurelist_clustered.csv'

### End user input

samplelist=pd.read_csv(data_dir+sample_list_name)
#samplefiles=samplelist.loc[samplelist['Sample type']=='sample','File']
#samplelist=samplelist[samplelist['Sample type']=='sample']

#Load MS data from sample list as MSfiles dictionary (keys=file name, values= parser objects)
MSfiles={}
outputs=[]
for file in samplelist['File']:
    print(file)
    parser = rawFileReader.ImportMassSpectraThermoMSFileReader(data_dir+file)
    MSfiles[file]=parser
    outputs.append(pd.read_csv(data_dir+file.replace('.raw','_assigned.csv')))

allresults=pd.concat(outputs,ignore_index=True)


#allresults['File']=allresults['file']
#allresults=allresults[allresults['m/z']<800]
#allresults=allresults[allresults['S/N']>3]
allresults=allresults[allresults['Time']>1]
#allresults=allresults[allresults['Time']<35]
allresults.to_csv(data_dir+'BATS_allresults.csv')

allresults.loc[allresults['Is Isotopologue']==1,'Molecular class']='Isotope'
allresults.loc[allresults['Heteroatom Class']=='unassigned','Molecular class']='unassigned'
allresults.loc[allresults['Molecular class']=='CH','Molecular class']='CHO'
allresults.loc[allresults['Molecular class']=='CHN','Molecular class']='CHON'


# Calculate atomic stoichiometries and Nominal Oxidation State of Carbon (NOSC)
results=allresults[allresults['Is Isotopologue']==0].fillna(0)
#results=results[(results['Molecular class'].isin(['CHO','CHON']))]

#results['O/C']=results['O']/results['C']
#results['H/C']=results['H']/results['C']
results['N/C']=results['N']/results['C']
results['NOSC'] =  4 -(4*results['C'] + results['H'] - 3*results['N'] - 2*results['O'])/results['C']


#### Plot library assignments over time

assign_summary=[]
for time in allresults['Time'].unique():
    current={}
    current['Time']=time
    for mol_class in allresults['Molecular class'].unique():
        current[mol_class]=len(allresults[(allresults['Molecular class']==mol_class) & (allresults['Time']==time)])
    assign_summary.append(current)

df=pd.DataFrame(assign_summary)

#Create a list of all unique features and their abundances across samples.
uniquelist=[]
for time in results.Time.unique():
    current=results[results.Time==time]
    current=current.sort_values(by=['Peak Height'],ascending=False)
    currentunique=current.drop_duplicates(subset=['Molecular Formula'])
    currentunique=currentunique.set_index(['Molecular Formula'],drop=False)
    for file in results['File'].unique():
        current_file=current[current['File']==file]
        current_file=current_file.rename(columns={'Peak Height':file})
        current_file=current_file.set_index(['Molecular Formula'],drop=False)
        currentunique=currentunique.join(current_file[file])
    uniquelist.append(currentunique)

featurelist=pd.concat(uniquelist,ignore_index=True)

print('All assignments:', len(results))
print("All Unique results: " + str(len(featurelist)))


### Performs statistical tests to evaluate reference classes
featurelist['rank class']=featurelist['Molecular Formula'].str.replace(' ', '',regex=True).str.replace('C','',regex=True).str.replace('H','',regex=True).str.replace('O','',regex=True).str.replace('\d+', '',1,regex=True)
featurelist['class flag']=0

ref_class='CHO'
filterby='rank class'

for c in featurelist[filterby].unique():
    #Calculate adjusted p-value from rank sum test (probability that distributions have same mean)
    ranksum=ranksums(featurelist[featurelist[filterby]==c]['m/z Error (ppm)'],
                   featurelist[featurelist['Molecular class']==ref_class]['m/z Error (ppm)'])
    class_size=len(featurelist[featurelist[filterby]==c])

    #Calculate statistic from (how similar distributions are)
    ks=ks_2samp(featurelist[featurelist[filterby]==c]['m/z Error (ppm)'],
                   featurelist[featurelist['Molecular class']==ref_class]['m/z Error (ppm)'])

    if class_size>9:
        if ranksum.pvalue*len(featurelist['rank class'].unique())>0.001:
            if ks.statistic<0.5:
                featurelist.loc[featurelist[filterby]==c,'class flag']=1


results1=featurelist[featurelist['class flag']==1]
results2=featurelist[featurelist['class flag']==0]

featurelist.loc[featurelist['class flag']==0,'Molecular class']=='Unassigned'

print(len(results1))
print(len(results2))

featurelist=featurelist[featurelist['class flag']==1]

#### Plot and save error distribution figure

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
fig.set_size_inches(8, 5)

sns.scatterplot(x='m/z',y='m/z Error (ppm)',hue='Molecular class',data=results1,ax=ax1, edgecolor='none')
ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
ax1.set_title('a', fontweight='bold', loc='left')

sns.kdeplot(x='m/z Error (ppm)',data=results1,hue='Molecular class',ax=ax2,legend=False,common_norm=False)
ax2.set_title('b', fontweight='bold', loc='left')

sns.scatterplot(x='m/z',y='m/z Error (ppm)',hue='Molecular class',data=results2,ax=ax3, edgecolor='none')
ax3.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
ax3.set_title('c', fontweight='bold', loc='left')

sns.kdeplot(x='m/z Error (ppm)',data=results2,hue='Molecular class',ax=ax4,legend=False)
ax4.set_title('d', fontweight='bold', loc='left')

fig.tight_layout()

fig.savefig(data_dir+'CoreLCMS_Fig_errorpot.pdf',dpi=300,format='pdf')


### (4) Feature filtering (abundance, blank subtraction, minimum samples) and metric plots

#remove low abundance hits
featurelist=featurelist[featurelist["S/N"]>3]
print("Unique results, S/N>3: " + str(len(featurelist)))

#remove features detected in the blank within 50% of the max intensity. 
featurelist['blank']=featurelist['RMB_190828_BATS24_blnk.raw'].fillna(0)/featurelist['Peak Height']
featurelist=featurelist[featurelist['blank']<0.5]
print("Unique results, blank subtracted: " + str(len(featurelist)))

#remove hits that don't appear in a minimum of 5 samples:
featurelist['occurrence']=featurelist[results['File'].unique()].gt(0).sum(axis=1)
featurelist=featurelist[featurelist['occurrence']>4.5]

print("Unique results, min thresh: " + str(len(featurelist)))
print("Unique molecular formula: " + str(len(featurelist['Molecular Formula'].unique())))

# Dispersity index calculations based on extracted ion chromatograms of each feature


#Calculate Dispersity Index. 
EIC={}
for file in featurelist['File'].unique():
    masses=featurelist[featurelist['File']==file]['m/z'].unique().tolist()
    EIC[file]=MSfiles[file].get_eics(target_mzs=masses,tic_data={},peak_detection=False,smooth=False)


dispersity=[]
for ind in featurelist.index:
    current=featurelist.loc[ind]
    time=[0,2]+current.Time
    file=current.File
    mass=current['m/z']
    chroma=pd.DataFrame({'EIC':EIC[file][0][mass].eic,'time':EIC[file][0][mass].time})
    chroma=chroma[chroma['time'].between(time[0],time[1])]
    chroma=chroma.sort_values(by='EIC',ascending=False)
    d=chroma[chroma.cumsum()['EIC']<0.5*chroma.sum()['EIC']].time.std()
    dispersity.append(d)

featurelist['Dispersity']=dispersity


#Generate plot of m/z error and molecular stoichiometry for all features (Fig S4)

fig, axs = plt.subplot_mosaic([['a','b','c']], figsize=(8,2.5), constrained_layout=True)

sns.scatterplot(x='m/z',y='m/z Error (ppm)',hue='Molecular class',data=featurelist,ax=axs['a'],s=10, legend=False)
#axs['a'].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
axs['a'].set_title('a', fontweight='bold', loc='left')

sns.scatterplot(x='O/C',y='H/C',hue='Molecular class',data=featurelist,ax=axs['b'],s=10, legend=False)
#axs['a'].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
axs['b'].set_title('b', fontweight='bold', loc='left')

sns.scatterplot(x='m/z',y='DBE',hue='Molecular class',data=featurelist,ax=axs['c'],s=10)
axs['c'].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., frameon=False)
axs['c'].set_title('c', fontweight='bold', loc='left')


fig.savefig(data_dir+'CoreLCMS_FigS4.pdf',dpi=300,format='pdf')


#Generate bar plot of formula assigned peaks and molecular stoichiometries over time (Fig 2)
fig, axs = plt.subplot_mosaic([['a','b'],['a','c']], figsize=(8,4), constrained_layout=True)

df.plot.bar(x='Time',y=df.columns[1:],stacked=True,ylabel='Peaks',ax=axs['a'])
axs['a'].set_xticklabels(axs['a'].get_xticklabels(),rotation=0)
axs['a'].set_title('a', fontweight='bold', loc='left')
#axs['a'].set_ylim(0,30000)
axs['a'].set(xlabel='Time (min)')

sns.violinplot(x="Time", y="O/C", data=featurelist, ax=axs['b'], legend=False, color='skyblue')
axs['b'].set(xlabel=None)
axs['b'].tick_params(right=True)
axs['b'].set_title('b', fontweight='bold', loc='left')

sns.violinplot(x="Time", y="H/C", data=featurelist, ax=axs['c'], legend=False, color='skyblue')
axs['c'].set(xlabel='Time (min)')
axs['c'].tick_params(right=True)
#axs['c'].set_title('c', fontweight='bold', loc='left')
axs['c'].sharex(axs['b'])

fig.savefig(data_dir+'CoreLCMS_Fig2.pdf',dpi=300,format='pdf')

featurelist.to_csv(data_dir+featurelist_file)
### (6) Hierarchical Clustering
### Generate clusters of ms features across depth.
plt.show()