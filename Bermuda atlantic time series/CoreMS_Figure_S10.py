# Python script for generating Figure S10
# RMB update 6/02/2023
# Contributors: Yuri Corilo, Will Kew, Christian Dewey, Rene Boiteau

########### Import the os module
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

hue_order = ['1 SLDOM', '2 SRDOM', '3 RDOM']


# Generate figure
sns.set_palette("colorblind",5)

fig, axs = plt.subplot_mosaic([['a','b','c']], figsize=(8,4), constrained_layout=True)

sns.scatterplot(x='C',y='H', data=clusteredlist[clusteredlist['cluster']=='1 SLDOM'], ax=axs['a'],legend=False, s=8)
axs['a'].set_title('SLDOM', fontweight='bold', loc='left')
axs['a'].set(xlim=(0,50),ylim=(0,100))

sns.scatterplot(x='C',y='H', data=clusteredlist[clusteredlist['cluster']=='2 SRDOM'], ax=axs['b'],legend=False, s=8)
axs['b'].set_title('SRDOM', fontweight='bold', loc='left')
axs['b'].set(xlim=(0,50),ylim=(0,100))

sns.scatterplot(x='C',y='H', data=clusteredlist[clusteredlist['cluster']=='3 RDOM'], ax=axs['c'],legend=False, s=8)
axs['c'].set_title('RDOM', fontweight='bold', loc='left')
axs['c'].set(xlim=(0,50),ylim=(0,100))

# Generate figure
fig, axs = plt.subplot_mosaic([['a','b','c']], figsize=(8,4), constrained_layout=True)

sns.stripplot(x='C',y='H', data=clusteredlist[(clusteredlist['O']==7) & (clusteredlist['N']==0)], native_scale=True, hue='cluster', hue_order=hue_order, ax=axs['a'],legend=False, s=2)
axs['a'].set_title('O7', fontweight='bold', loc='left')
axs['a'].set(xlim=(0,50),ylim=(0,100))


sns.stripplot(x='C',y='H', data=clusteredlist[(clusteredlist['O']==7) & (clusteredlist['N']==1)], native_scale=True, hue='cluster', hue_order=hue_order, ax=axs['b'],legend=False, s=2)
axs['b'].set_title('O7 N1', fontweight='bold', loc='left')
axs['b'].set(xlim=(0,50),ylim=(0,100))

sns.stripplot(x='C',y='H', data=clusteredlist[(clusteredlist['O']==8) & (clusteredlist['N']==1)], native_scale=True, hue='cluster', hue_order=hue_order, ax=axs['c'],legend=True, s=2)
axs['c'].set_title('O8 N1', fontweight='bold', loc='left')
axs['c'].set(xlim=(0,50),ylim=(0,100))


fig.savefig(data_dir+'CoreLCMS_Fig_CH.pdf',dpi=300,format='pdf')
