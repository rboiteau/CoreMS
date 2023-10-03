# Import the os module
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import KMeans
from scipy.stats import ranksums
from matplotlib.colors import LogNorm, Normalize


### Set data directory and file names here

data_dir = '/CoreMS/usrdata/'
sample_list_name='BATS_sample_list.csv' #Sample list must contain column with header 'File'
featurelist_file='BATS_featurelist.csv'
clustered_featurelist_file='BATS_featurelist_clustered.csv'

n_clust=4

### End user input

featurelist=pd.read_csv(data_dir+featurelist_file)

samplelist=pd.read_csv(data_dir+sample_list_name)
samplelist['Depth']=samplelist['Depth (m)']

depth=[]
for file in featurelist['File'].unique():
    d=samplelist[samplelist['File']==file]['Depth'].iloc[0]
    depth.append(d)


### (6) Hierarchical Clustering
### Generate clusters of ms features across depth.


norm_abundances=featurelist[featurelist['File'].unique()].fillna(0)
norm_abundances=norm_abundances.div(norm_abundances.max(axis=1),axis=0)

### Make a heat map w/ heirarchical clustering
#sns.clustermap(norm_abundances,row_cluster=True,cmap='mako',method='average')


Kmean = KMeans(n_clusters=n_clust, n_init='auto',random_state=1)
Kmean.fit(norm_abundances)

featurelist['cluster']=Kmean.labels_


featurelist.to_csv(data_dir+clustered_featurelist_file)
