# Import the os module
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering
from scipy.stats import ranksums
from matplotlib.colors import LogNorm, Normalize

file_location='/Users/boiteaur/Desktop/Major projects/Bermuda Atlantic Time Series data processing/Thermo RAW data/'
cluster_results_file='BATS_clustered_results_kmeans.csv'

results_clustered=pd.read_csv(file_location+cluster_results_file)
results_clustered=results_clustered.fillna(0)

sample_list_name='BATS_sample_list.csv' #Sample list must contain column with header 'File'
samplelist=pd.read_csv(file_location+sample_list_name)
samplelist['Depth']=samplelist['Depth (m)']

clustermethod='average'

depth=[]
for file in results_clustered['File'].unique():
    d=samplelist[samplelist['File']==file]['Depth'].iloc[0]
    depth.append(d)


### (6) Hierarchical Clustering
### Module defines clusters based on depth distribution and plots density distribution of properties of each cluster:

#Define clusters
#results_clustered.cluster[results_clustered.cluster==0]='out'
#results_clustered.cluster[results_clustered.cluster==1]='3 RDOM'
#results_clustered.cluster[results_clustered.cluster==2]='1 SLDOM'
#results_clustered.cluster[results_clustered.cluster==3]='2 SRDOM'


results_clustered.cluster[results_clustered.cluster==0]='1 SLDOM'
results_clustered.cluster[results_clustered.cluster==1]='2 SRDOM'
results_clustered.cluster[results_clustered.cluster==2]='3 RDOM'
results_clustered.cluster[results_clustered.cluster==3]='out'


# Save clustered results (Table S1)
results_clustered.to_csv(file_location+'TableS1_rev_clustered_results.csv')

#print(len(results_clustered[results_clustered.cluster=='out']))

# Discard out cluster
results_clustered=results_clustered[results_clustered.cluster!='out']
results_clustered=results_clustered.sort_values(by='cluster')

# Generate figures (Fig. 4)
sns.set_palette("colorblind",5)

# Generate violin plots version
fig2, axs = plt.subplot_mosaic([['a','b','c'],['d','e','f']], figsize=(8,5), constrained_layout=True)
sns.violinplot(x='H/C',y='cluster', data=results_clustered, ax=axs['a'], common_norm=False,legend=False)
axs['a'].set_title('a', fontweight='bold', loc='left')
sns.violinplot(x='O/C',y='cluster', data=results_clustered, ax=axs['b'], common_norm=False,legend=False)
axs['b'].set_title('b', fontweight='bold', loc='left')
sns.violinplot(x='N/C',y='cluster', data=results_clustered, ax=axs['c'], common_norm=False,legend=False)
axs['c'].set_title('c', fontweight='bold', loc='left')
axs['c'].set_xlim(0,0.3)
sns.violinplot(x='m/z',y='cluster', data=results_clustered, ax=axs['d'], common_norm=False,legend=False)
axs['d'].set_title('d', fontweight='bold', loc='left')
axs['d'].set_xlabel("$\it{m/z}$")
sns.violinplot(x='NOSC',y='cluster', data=results_clustered, ax=axs['e'], common_norm=False,legend=False)
axs['e'].set_title('e', fontweight='bold', loc='left')
sns.violinplot(x='Time',y='cluster', data=results_clustered, ax=axs['f'], common_norm=False,legend=False)
axs['f'].set_title('f', fontweight='bold', loc='left')

# Generate violin plots version
fig2, axs = plt.subplot_mosaic([['a','b','c'],['d','e','f']], figsize=(8,5), constrained_layout=True)
sns.boxplot(x='H/C',y='cluster', data=results_clustered, ax=axs['a'])
axs['a'].set_title('a', fontweight='bold', loc='left')
sns.boxplot(x='O/C',y='cluster', data=results_clustered, ax=axs['b'])
axs['b'].set_title('b', fontweight='bold', loc='left')
sns.boxplot(x='N/C',y='cluster', data=results_clustered, ax=axs['c'])
axs['c'].set_title('c', fontweight='bold', loc='left')
axs['c'].set_xlim(0,0.3)
sns.boxplot(x='m/z',y='cluster', data=results_clustered, ax=axs['d'])
axs['d'].set_title('d', fontweight='bold', loc='left')
axs['d'].set_xlabel("$\it{m/z}$")
sns.boxplot(x='NOSC',y='cluster', data=results_clustered, ax=axs['e'])
axs['e'].set_title('e', fontweight='bold', loc='left')
sns.boxplot(x='Time',y='cluster', data=results_clustered, ax=axs['f'])
axs['f'].set_title('f', fontweight='bold', loc='left')


### (7) Cluster analysis
### Module determines and plots the mean depth distribution of each cluster

norm_abundances=results_clustered[results_clustered['File'].unique()].fillna(0)
norm_abundances=norm_abundances.div(norm_abundances.max(axis=1),axis=0)

depth=[]
for file in results_clustered['File'].unique():
    d=samplelist[samplelist['File']==file]['Depth'].iloc[0]
    depth.append(d)

clustered_results=[]
for i in results_clustered['cluster'].unique():
    current=pd.DataFrame({'abundance':norm_abundances[results_clustered['cluster']==i].sum(axis=0)/len(results_clustered[results_clustered['cluster']==i]),'Depth':depth,'cluster':i})
    clustered_results.append(current)

clustered_results=pd.concat(clustered_results)

h = sns.relplot(data=clustered_results, col='cluster', x='abundance', y='Depth', kind='scatter', height=3)
for ax in h.axes[0]:
    ax.set_ylim(1000,0)
    ax.set_xlabel('Normalized Abundance')

h.savefig(file_location+'CoreLCMS_FigS6.eps',dpi=300,format='eps')
h.savefig(file_location+'CoreLCMS_FigS6.pdf',dpi=300,format='pdf')

#Print number of features in each cluster
print(len(results_clustered))

#Create and saves heatmap for each cluster (Fig. 5)

for cluster in results_clustered['cluster'].unique():

    current=results_clustered[results_clustered['cluster']==cluster]
    current=current[results_clustered['File'].unique()].fillna(0)
    current=current.div(current.max(axis=1),axis=0)
    print(cluster)
    print(len(current))
    clusterplot=norm_abundances
    current=current.transpose()
    current['Depth']=depth
    current=current.sort_values(by='Depth')
    current=current.drop(['Depth'],axis=1)
    h=sns.clustermap(current,row_cluster=False,cmap='mako',method=clustermethod)
    h.savefig(file_location+'CoreLCMS_Fig3'+cluster+'.eps',dpi=300,format='eps')


plt.show()