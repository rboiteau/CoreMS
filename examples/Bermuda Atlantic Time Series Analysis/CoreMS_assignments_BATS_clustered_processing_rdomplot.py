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

clustermethod='ward'

depth=[]
for file in results_clustered['File'].unique():
    d=samplelist[samplelist['File']==file]['Depth'].iloc[0]
    depth.append(d)


### (6) Hierarchical Clustering
### Module defines clusters based on depth distribution and plots density distribution of properties of each cluster:

#Define clusters
results_clustered.cluster[results_clustered.cluster==0]='out'
results_clustered.cluster[results_clustered.cluster==1]='RDOM'
results_clustered.cluster[results_clustered.cluster==2]='1 SLDOM'
results_clustered.cluster[results_clustered.cluster==3]='2 SRDOM'


# Save clustered results (Table S1)
#results_clustered.to_csv(file_location+'TableS1_rev_clustered_results.csv')

#print(len(results_clustered[results_clustered.cluster=='out']))

# Discard out cluster
#results_clustered=results_clustered[results_clustered.cluster=='RDOM']
results_clustered=results_clustered[results_clustered.cluster!='out']

results_clustered=results_clustered[results_clustered['Molecular class'].isin(['CHO','CHON'])]

#Calculate kendrick mass defect.
results_clustered['KM']=results_clustered['m/z']*14/14.01565
results_clustered['KMD']=results_clustered['KM'].round(0)-results_clustered['KM']

# Generate figures (Fig. 4)
sns.set_palette("colorblind",5)

fig2, axs = plt.subplot_mosaic([['a','b']], figsize=(6,3), constrained_layout=True)
sns.scatterplot(x='O/C',y='H/C',hue='Molecular class', data=results_clustered, ax=axs['a'], legend=True)
axs['a'].set_title('a', fontweight='bold', loc='left')
sns.scatterplot(x='KM',y='KMD',hue='Molecular class', data=results_clustered, ax=axs['b'], legend=False)
axs['b'].set_title('b', fontweight='bold', loc='left')


plt.show()