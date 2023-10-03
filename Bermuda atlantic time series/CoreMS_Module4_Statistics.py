# Import the os module
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering
from scipy.stats import ranksums
from matplotlib.colors import LogNorm, Normalize


### Set data directory and file names here

data_dir = '/CoreMS/usrdata/'
sample_list_name='BATS_sample_list.csv' #Sample list must contain column with header 'File'
featurelist_file='BATS_featurelist.csv'
clustered_featurelist_file='BATS_featurelist_clustered.csv'

### End user input


clusteredlist=pd.read_csv(data_dir+clustered_featurelist_file).fillna(0)
samplelist=pd.read_csv(data_dir+sample_list_name)
samplelist['Depth']=samplelist['Depth (m)']
depth=[]
for file in clusteredlist['File'].unique():
    d=samplelist[samplelist['File']==file]['Depth'].iloc[0]
    depth.append(d)


clusteredlist['NOSC'] =  4 -(4*clusteredlist['C'] + clusteredlist['H'] - 3*clusteredlist['N'] - 2*clusteredlist['O'])/clusteredlist['C']

### (6) Hierarchical Clustering
### Module defines clusters based on depth distribution and plots density distribution of properties of each cluster:

#Define clusters
clusteredlist.cluster[clusteredlist.cluster==0]='out'
clusteredlist.cluster[clusteredlist.cluster==1]='1 SLDOM'
clusteredlist.cluster[clusteredlist.cluster==2]='2 SRDOM'
clusteredlist.cluster[clusteredlist.cluster==3]='3 RDOM'


# Save clustered results (Table S1)
clusteredlist.to_csv(data_dir+'TableS1_featurelist.csv')

print(len(clusteredlist[clusteredlist.cluster=='out']))

# Discard out cluster
clusteredlist=clusteredlist[clusteredlist.cluster!='out']
clusteredlist=clusteredlist.sort_values(by='cluster')

# Generate figures (Fig. 4)
sns.set_palette("colorblind",5)

# Generate violin plots version
fig, axs = plt.subplot_mosaic([['a','b','c'],['d','e','f']], figsize=(8,5), constrained_layout=True)
sns.violinplot(x='H/C',y='cluster', data=clusteredlist, ax=axs['a'], common_norm=False,legend=False)
axs['a'].set_title('a', fontweight='bold', loc='left')
sns.violinplot(x='O/C',y='cluster', data=clusteredlist, ax=axs['b'], common_norm=False,legend=False)
axs['b'].set_title('b', fontweight='bold', loc='left')
axs['b'].set(ylabel=None,yticklabels=[])
sns.violinplot(x='NOSC',y='cluster', data=clusteredlist, ax=axs['c'], common_norm=False,legend=False)
axs['c'].set_title('c', fontweight='bold', loc='left')
axs['c'].set(ylabel=None,yticklabels=[])
sns.violinplot(x='m/z',y='cluster', data=clusteredlist, ax=axs['d'], common_norm=False,legend=False)
axs['d'].set_title('d', fontweight='bold', loc='left')
axs['d'].set_xlabel("$\it{m/z}$")
sns.violinplot(x='Dispersity',y='cluster', data=clusteredlist, ax=axs['e'], common_norm=False,legend=False)
axs['e'].set_title('e', fontweight='bold', loc='left')
axs['e'].set(ylabel=None,yticklabels=[])
sns.violinplot(x='Time',y='cluster', data=clusteredlist, ax=axs['f'], common_norm=False,legend=False)
axs['f'].set_title('f', fontweight='bold', loc='left')
axs['f'].set(ylabel=None,yticklabels=[],xlabel='Retention time (min)')

fig.savefig(data_dir+'CoreLCMS_Fig4.pdf',dpi=300,format='pdf')



### Performs statistical tests to determine significant compositional differences between across clusters. 


stat='H/C'
g1='3 RDOM'
g2='1 SLDOM'
print(stat)
print(g1)
print(clusteredlist[clusteredlist['cluster']==g1][stat].mean())
print(clusteredlist[clusteredlist['cluster']==g2][stat].mean())
print(ranksums(clusteredlist[clusteredlist['cluster']=='3 RDOM'][stat],clusteredlist[clusteredlist['cluster']=='1 SLDOM'][stat]))

stat='O/C'
g1='3 RDOM'
g2='1 SLDOM'
print(stat)
print(g1)
print(clusteredlist[clusteredlist['cluster']==g1][stat].mean())
print(clusteredlist[clusteredlist['cluster']==g2][stat].mean())
print(ranksums(clusteredlist[clusteredlist['cluster']=='3 RDOM'][stat],clusteredlist[clusteredlist['cluster']=='1 SLDOM'][stat]))

stat='N/C'
g1='3 RDOM'
g2='1 SLDOM'
print(stat)
print(g1)
print(clusteredlist[clusteredlist['cluster']==g1][stat].mean())
print(clusteredlist[clusteredlist['cluster']==g2][stat].mean())
print(ranksums(clusteredlist[clusteredlist['cluster']=='3 RDOM'][stat],clusteredlist[clusteredlist['cluster']=='1 SLDOM'][stat]))

stat='NOSC'
g1='3 RDOM'
g2='1 SLDOM'
print(stat)
print(g1)
print(clusteredlist[clusteredlist['cluster']==g1][stat].mean())
print(clusteredlist[clusteredlist['cluster']==g2][stat].mean())
print(ranksums(clusteredlist[clusteredlist['cluster']=='3 RDOM'][stat],clusteredlist[clusteredlist['cluster']=='1 SLDOM'][stat]))

stat='m/z'
g1='3 RDOM'
g2='1 SLDOM'
print(stat)
print(g1)
print(clusteredlist[clusteredlist['cluster']==g1][stat].mean())
print(clusteredlist[clusteredlist['cluster']==g2][stat].mean())
print(ranksums(clusteredlist[clusteredlist['cluster']=='3 RDOM'][stat],clusteredlist[clusteredlist['cluster']=='1 SLDOM'][stat]))

stat='Dispersity'
g1='3 RDOM'
g2='1 SLDOM'
print(stat)
print(g1)
print(clusteredlist[clusteredlist['cluster']==g1][stat].mean())
print(clusteredlist[clusteredlist['cluster']==g2][stat].mean())
print(ranksums(clusteredlist[clusteredlist['cluster']=='3 RDOM'][stat],clusteredlist[clusteredlist['cluster']=='1 SLDOM'][stat]))

