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

clusteredlist=clusteredlist[clusteredlist['S/N']>3]

### (6) Hierarchical Clustering
### Module defines clusters based on depth distribution and plots density distribution of properties of each cluster:

#Define clusters
clusteredlist.cluster[clusteredlist.cluster==0]='out'
clusteredlist.cluster[clusteredlist.cluster==1]='1 SLDOM'
clusteredlist.cluster[clusteredlist.cluster==2]='2 SRDOM'
clusteredlist.cluster[clusteredlist.cluster==3]='3 RDOM'


# Save clustered results (Table S1)
clusteredlist.to_csv(data_dir+'TableS1_featurelist.csv')


# Discard out cluster
clusteredlist=clusteredlist[clusteredlist.cluster!='out']
#clusteredlist=clusteredlist.sort_values(by='cluster')
clusteredlist=clusteredlist.sort_values(by='N')
clusteredlist=clusteredlist.sort_values(by='O')
#clusteredlist=clusteredlist.sort_values(by='cluster')





print(len(clusteredlist))

# Generate figures (Fig. 4)
sns.set_palette("colorblind",5)

# Generate violin plots version
fig, axs = plt.subplot_mosaic([['a','b','c']], figsize=(8,4), constrained_layout=True)
#fig, axs = plt.subplot_mosaic([['a']], figsize=(8,5), constrained_layout=True)

#sns.scatterplot(x='m/z',y='DBE', data=clusteredlist, ax=axs['a'], hue='cluster',legend=False, s=8)
#axs['a'].set_title('a', fontweight='bold', loc='left')
#sns.scatterplot(x='m/z',y='H/C', data=clusteredlist, ax=axs['b'], hue='cluster',legend=True, s=8)
#axs['b'].set_title('c', fontweight='bold', loc='left')

#fig=sns.jointplot(x='m/z',y='H/C', data=clusteredlist, hue='cluster',legend=True, s=8)

hue_order = ['1 SLDOM', '2 SRDOM', '3 RDOM']

sns.histplot(x='Heteroatom Class', data=clusteredlist[clusteredlist['N']==0], ax=axs['a'], hue='cluster', hue_order=hue_order, legend=False,multiple='dodge')
axs['a'].tick_params(axis='x', rotation=90)
axs['a'].set(xlabel=' ')
axs['a'].set(ylim=(0,500))

sns.histplot(x='Heteroatom Class', data=clusteredlist[clusteredlist['N']==1], ax=axs['b'], hue='cluster', hue_order=hue_order, legend=False,multiple='dodge')
axs['b'].tick_params(axis='x', rotation=90)
axs['b'].set(ylabel=None,yticklabels=[],xlabel='Heteroatom class')
axs['b'].set(ylim=(0,500))

sns.histplot(x='Heteroatom Class', data=clusteredlist[clusteredlist['N']==2], ax=axs['c'], hue='cluster', hue_order=hue_order, legend=True,multiple='dodge')
axs['c'].tick_params(axis='x', rotation=90)
axs['c'].set(ylabel=None,yticklabels=[],xlabel=' ')
axs['c'].set(ylim=(0,500))


fig.savefig(data_dir+'CoreLCMS_Fig_classes.pdf',dpi=300,format='pdf')
