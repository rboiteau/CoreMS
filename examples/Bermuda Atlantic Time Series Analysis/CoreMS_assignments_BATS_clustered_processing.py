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

results_clustered=pd.read_csv(file_location+'BATS_clustered_results.csv')

sample_list_name='BATS_sample_list.csv' #Sample list must contain column with header 'File'
samplelist=pd.read_csv(file_location+sample_list_name)
samplelist['Depth']=samplelist['Depth (m)']

depth=[]
for file in results_clustered['File'].unique():
    d=samplelist[samplelist['File']==file]['Depth'].iloc[0]
    depth.append(d)


### (6) Hierarchical Clustering
### Module defines clusters based on depth distribution and plots density distribution of properties of each cluster:

#Define clusters
results_clustered.cluster[results_clustered.cluster==0]='3.SRDOM'
results_clustered.cluster[results_clustered.cluster==1]='2.SLDOM'
results_clustered.cluster[results_clustered.cluster==2]='out'
results_clustered.cluster[results_clustered.cluster==3]='out'
results_clustered.cluster[results_clustered.cluster==4]='4.RDOM'
results_clustered.cluster[results_clustered.cluster==5]='1.SLDOM'
results_clustered.cluster[results_clustered.cluster==6]='4.RDOM'
results_clustered.cluster[results_clustered.cluster==7]='3.SRDOM'
results_clustered.cluster[results_clustered.cluster==8]='4.RDOM'

# Save clustered results (Table S1)
results_clustered.to_csv(file_location+'TableS1_rev_clustered_results.csv')

# Discard out cluster
results_clustered=results_clustered[results_clustered.cluster!='out']
results_clustered=results_clustered.sort_values(by='cluster')

# Generate figures (Fig. 4)
sns.set_palette("colorblind",5)

fig2, axs = plt.subplot_mosaic([['a','b','c'],['d','e','f']], figsize=(8,5), constrained_layout=True)
sns.kdeplot(x='H/C',hue='cluster', data=results_clustered, ax=axs['a'], common_norm=False,legend=False)
axs['a'].set_title('a', fontweight='bold', loc='left')
sns.kdeplot(x='O/C',hue='cluster', data=results_clustered, ax=axs['b'], common_norm=False,legend=False)
axs['b'].set_title('b', fontweight='bold', loc='left')
sns.kdeplot(x='N/C',hue='cluster', data=results_clustered, ax=axs['c'], common_norm=False,legend=True)
axs['c'].set_title('c', fontweight='bold', loc='left')
sns.kdeplot(x='m/z',hue='cluster', data=results_clustered, ax=axs['d'], common_norm=False,legend=False)
axs['d'].set_title('d', fontweight='bold', loc='left')
axs['d'].set_xlabel("$\it{m/z}$")
sns.kdeplot(x='NOSC',hue='cluster', data=results_clustered, ax=axs['e'], common_norm=False,legend=False)
axs['e'].set_title('e', fontweight='bold', loc='left')
sns.kdeplot(x='Time',hue='cluster', data=results_clustered, ax=axs['f'], common_norm=False,legend=False)
axs['f'].set_title('f', fontweight='bold', loc='left')

fig2.savefig(file_location+'CoreLCMS_Fig4.eps',dpi=300,format='eps')
fig2.savefig(file_location+'CoreLCMS_Fig4.pdf',dpi=300,format='pdf')

# Generate violin plots version
fig2, axs = plt.subplot_mosaic([['a','b','c'],['d','e','f']], figsize=(8,5), constrained_layout=True)
sns.violinplot(x='H/C',y='cluster', data=results_clustered, ax=axs['a'], common_norm=False,legend=False)
axs['a'].set_title('a', fontweight='bold', loc='left')
sns.violinplot(x='O/C',y='cluster', data=results_clustered, ax=axs['b'], common_norm=False,legend=False)
axs['b'].set_title('b', fontweight='bold', loc='left')
sns.violinplot(x='N/C',y='cluster', data=results_clustered, ax=axs['c'], common_norm=False,legend=False)
axs['c'].set_title('c', fontweight='bold', loc='left')
sns.violinplot(x='m/z',y='cluster', data=results_clustered, ax=axs['d'], common_norm=False,legend=False)
axs['d'].set_title('d', fontweight='bold', loc='left')
axs['d'].set_xlabel("$\it{m/z}$")
sns.violinplot(x='NOSC',y='cluster', data=results_clustered, ax=axs['e'], common_norm=False,legend=False)
axs['e'].set_title('e', fontweight='bold', loc='left')
sns.violinplot(x='Time',y='cluster', data=results_clustered, ax=axs['f'], common_norm=False,legend=False)
axs['f'].set_title('f', fontweight='bold', loc='left')

plt.show()