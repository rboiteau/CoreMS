# Import the os module
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

file_location='/Users/boiteaur/Desktop/Major projects/Bermuda Atlantic Time Series data processing/Thermo RAW data/'

allresults=pd.read_csv(file_location+'BATS_allfiles_assigned_results.csv')

#sample_list_name='BATS_sample_list.csv' #Sample list must contain column with header 'File'

#Load MS data from sample list as MSfiles dictionary (keys=file name, values= parser objects)
#samplelist=pd.read_csv(file_location+sample_list_name)
#MSfiles={}
#for file in samplelist['File']:
#    parser = rawFileReader.ImportMassSpectraThermoMSFileReader(file_location+file)
#    MSfiles[file]=parser
allresults=allresults[allresults['Time']<26]
allresults=allresults[allresults['Time']>3]

allresults['Molecular class'][allresults['N']>4]=='CHON>4'

results=allresults[allresults['Is Isotopologue']==0]
results=results[(results['Molecular class'].isin(['CHO','CHON']))]

# Calculate atomic stoichiometries and Nominal Oxidation State of Carbon (NOSC)
results['O/C']=results['O']/results['C']
results['H/C']=results['H']/results['C']
results['N/C']=results['N']/results['C']
results['NOSC'] =  4 -(4*results['C'] + results['H'] - 3*results['N'] - 2*results['O'])/results['C']

print('All assignments:', len(results))

'''
fig, ((ax1, ax2)) = plt.subplots(1,2)
fig.set_size_inches(10, 7)

#### Plot and save error distribution figure
sns.scatterplot(x='m/z',y='m/z Error (ppm)',hue='Molecular class',data=results,ax=ax1, edgecolor='none')
ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
ax1.set_xlim(200,800)
sns.kdeplot(x='m/z Error (ppm)',data=results,hue='Molecular class',ax=ax2,legend=False)
ax2.set_xlim(-0.3,0.3)
fig.tight_layout()
'''

#### Plot library assignments over time

assign_summary=[]
for time in allresults['Time'].unique():
    current={}
    current['Time']=time
    for mol_class in allresults['Molecular class'].unique():
        current[mol_class]=len(allresults[(allresults['Molecular class']==mol_class) & (allresults['Time']==time)])
    assign_summary.append(current)
    #mzdiff=result['m/z'].sort_values(ascending=True).diff().iloc[1:]/result['m/z'].sort_values(ascending=True).iloc[1:]*1E6

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

uniqueresults=pd.concat(uniquelist,ignore_index=True)
print("All Unique results: " + str(len(uniqueresults)))

### (4) Feature filtering (abundance, blank subtraction, minimum samples) and metric plots

#remove low abundance hits
uniqueresults=uniqueresults[uniqueresults["S/N"]>3]
print("Unique results, S/N>3: " + str(len(uniqueresults)))

#remove features detected in the blank within 50% of the max intensity. 
uniqueresults['blank']=uniqueresults['RMB_190828_BATS24_blnk.raw'].fillna(0)/uniqueresults['Peak Height']
uniqueresults=uniqueresults[uniqueresults['blank']<0.5]
print("Unique results, blank subtracted: " + str(len(uniqueresults)))

#remove hits that don't appear in a minimum of 5 samples:
uniqueresults['occurrence']=uniqueresults[results['File'].unique()].gt(0).sum(axis=1)
uniqueresults=uniqueresults[uniqueresults['occurrence']>4.5]

print("Unique results, min thresh: " + str(len(uniqueresults)))
print("Unique molecular formula: " + str(len(uniqueresults['Molecular Formula'].unique())))

#Generate plot of m/z error and molecular stoichiometry for all features (Fig S4)

fig, axs = plt.subplot_mosaic([['a','b','c']], figsize=(7.5,2.5), constrained_layout=True)

sns.scatterplot(x='m/z',y='m/z Error (ppm)',hue='Molecular class',data=uniqueresults,ax=axs['a'],s=10, legend=False)
#axs['a'].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
axs['a'].set_title('a', fontweight='bold', loc='left')

sns.scatterplot(x='O/C',y='H/C',hue='Molecular class',data=uniqueresults,ax=axs['b'],s=10, legend=False)
#axs['a'].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
axs['b'].set_title('b', fontweight='bold', loc='left')

sns.scatterplot(x='m/z',y='DBE',hue='Molecular class',data=uniqueresults,ax=axs['c'],s=10)
axs['c'].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., frameon=False)
axs['c'].set_title('c', fontweight='bold', loc='left')


fig.savefig(file_location+'CoreLCMS_FigS4.eps',dpi=300,format='eps')
fig.savefig(file_location+'CoreLCMS_FigS4.pdf',dpi=300,format='pdf')


#Generate bar plot of formula assigned peaks and molecular stoichiometries over time (Fig 2)
fig, axs = plt.subplot_mosaic([['a','b'],['a','c']], figsize=(7,4), constrained_layout=True)

df.plot.bar(x='Time',y=df.columns[1:],stacked=True,ylabel='Peaks',ax=axs['a'])
axs['a'].set_xticklabels(axs['a'].get_xticklabels(),rotation=0)
axs['a'].set_title('a', fontweight='bold', loc='left')
#axs['a'].set_ylim(0,30000)
axs['a'].set(xlabel='Time (min)')

sns.violinplot(x="Time", y="O/C", data=uniqueresults, ax=axs['b'], legend=False, color='skyblue')
axs['b'].set(xlabel=None)
axs['b'].tick_params(right=True)
axs['b'].set_title('b', fontweight='bold', loc='left')

sns.violinplot(x="Time", y="H/C", data=uniqueresults, ax=axs['c'], legend=False, color='skyblue')
axs['c'].set(xlabel='Time (min)')
axs['c'].tick_params(right=True)
#axs['c'].set_title('c', fontweight='bold', loc='left')
axs['c'].sharex(axs['b'])

fig.savefig(file_location+'CoreLCMS_Fig2.eps',dpi=300,format='eps')
fig.savefig(file_location+'CoreLCMS_Fig2.pdf',dpi=300,format='pdf')

uniqueresults.to_csv(file_location+'BATS_unique_results.csv')
### (6) Hierarchical Clustering
### Generate clusters of ms features across depth.
