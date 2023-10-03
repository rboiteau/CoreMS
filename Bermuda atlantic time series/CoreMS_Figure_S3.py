# Python script for generating Figure S3: Evaluation of molecular formula stoichiometry search parameters by mass error. 
# RMB update 6/02/2023

# Import the os module
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

file_location='/Users/boiteaur/Desktop/Major projects/Bermuda Atlantic Time Series data processing/Thermo RAW data/'

allresults=pd.read_csv(file_location+'BATSpooled_assigned_results.csv')
#allresults=allresults[allresults['Time']<28]
#allresults=allresults[allresults['Time']>2]

print(len(allresults))

allresults['Molecular class'][allresults['N']>5]='CHON>4'
allresults['Molecular class'][allresults['S']>1]='CHOS>2'

results=allresults[allresults['Is Isotopologue']==0]
#results=allresults[allresults['Molecular class']!='unassigned']

results_1=results[(results['Molecular class'].isin(['CHO','CHON']))]
results_2=results[(~results['Molecular class'].isin(['CHO','CHON']))]

print('All assignments:', len(results))
print('CHON assignments:', len(results_1))
print('Other assignments:', len(results_2))

#Generate figure
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
fig.set_size_inches(8, 5)

sns.scatterplot(x='m/z',y='m/z Error (ppm)',hue='Molecular class',data=results_1,ax=ax1, edgecolor='none')
ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
ax1.set_xlim(200,900)
sns.kdeplot(x='m/z Error (ppm)',data=results_1,hue='Molecular class',ax=ax2,legend=False)
ax2.set_xlim(-0.3,0.3)
fig.tight_layout()


sns.scatterplot(x='m/z',y='m/z Error (ppm)',hue='Molecular class',data=results_2,ax=ax3, edgecolor='none')
ax3.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
ax3.set_xlim(200,900)
sns.kdeplot(x='m/z Error (ppm)',data=results_2,hue='Molecular class',ax=ax4,legend=False)
ax4.set_xlim(-0.3,0.3)
fig.tight_layout()

fig.savefig(file_location+'CoreLCMS_FigS3.eps',dpi=300,format='eps')
fig.savefig(file_location+'CoreLCMS_FigS3.pdf',dpi=300,format='pdf')

plt.show()
