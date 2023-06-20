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

uniqueresults=pd.read_csv(file_location+'BATS_unique_results_round3.csv')

sample_list_name='BATS_sample_list.csv' #Sample list must contain column with header 'File'
samplelist=pd.read_csv(file_location+sample_list_name)
samplelist['Depth']=samplelist['Depth (m)']

depth=[]
for file in uniqueresults['File'].unique():
    d=samplelist[samplelist['File']==file]['Depth'].iloc[0]
    depth.append(d)


### (6) Hierarchical Clustering
### Generate clusters of ms features across depth.

#Cluster settings
clustermethod='ward'
nclusters=9

#Clustering functions

current=uniqueresults[uniqueresults['File'].unique()].fillna(0)
current=current.div(current.max(axis=1),axis=0)
#clusterplot=norm_abundances
current=current.transpose()
current['Depth']=depth
current=current.sort_values(by='Depth')
current=current.drop(['Depth'],axis=1)
sns.clustermap(current,row_cluster=False,cmap='mako',method=clustermethod)

plt.show()

###
results_clustered=uniqueresults

norm_abundances=results_clustered[uniqueresults['File'].unique()].fillna(0)
norm_abundances=norm_abundances.div(norm_abundances.max(axis=1),axis=0)


cluster = AgglomerativeClustering(n_clusters=nclusters,affinity='euclidean',linkage=clustermethod)
cluster.fit_predict(norm_abundances)

results_clustered['cluster']=cluster.labels_
print(len(results_clustered))

clustered_results=[]
for i in results_clustered['cluster'].unique():
    current=pd.DataFrame({'abundance':norm_abundances[results_clustered['cluster']==i].sum(axis=0)/len(results_clustered[results_clustered['cluster']==i]),'Depth':depth,'cluster':i})
    clustered_results.append(current)

clustered_results=pd.concat(clustered_results)
clustered_results=clustered_results.sort_values(by='Depth',ascending=True)

g = sns.relplot(data=clustered_results, col='cluster', x='abundance', y='Depth', kind='scatter')
for ax in g.axes[0]:
    ax.invert_yaxis()

plt.show()

results_clustered.to_csv(file_location+'BATS_clustered_results_round3.csv')
