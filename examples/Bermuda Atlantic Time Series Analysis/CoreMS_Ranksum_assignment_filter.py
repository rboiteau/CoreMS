#Statistical test to weed out bad assignment types. 
# Import other required modules
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from matplotlib.colors import LogNorm, Normalize
from scipy.stats import ranksums
from scipy.stats import kruskal
from scipy.stats import fligner
from scipy.stats import bartlett
from scipy.stats import levene
from scipy.stats import shapiro
from unidip import UniDip
from scipy.stats import ks_2samp

file_location='/Users/boiteaur/Desktop/Major projects/Bermuda Atlantic Time Series data processing/Thermo RAW data/'

allresults=pd.read_csv(file_location+'BATSpooled_assigned_results.csv')

print(len(allresults))
allresults=allresults[allresults['m/z']<900]
allresults=allresults[allresults['m/z Error (ppm)']<0.25]
allresults=allresults[allresults['m/z Error (ppm)']>-0.25]

#allresults=allresults[allresults['Time']>3]
allresults=allresults[allresults['S/N']>3]


### Performs statistical tests to determine significant compositional differences between across clusters. 
allresults['rank class']=allresults['Molecular Formula'].str.replace(' ', '',regex=True).str.replace('C','',regex=True).str.replace('H','',regex=True).str.replace('O','',regex=True).str.replace('\d+', '',1,regex=True)
allresults['class flag']=0


assignedresults=allresults[allresults['Is Isotopologue']==0]
#results=results[(results['Molecular class'].isin(['CHO','CHON']))]
ref_class='CHO'

filterby='rank class'
for c in assignedresults[filterby].unique():
    #Calculate adjusted p-value from rank sum test (probability that distributions have same mean)
    ranksum=ranksums(assignedresults[assignedresults[filterby]==c]['m/z Error (ppm)'],
                   assignedresults[assignedresults['Molecular class']==ref_class]['m/z Error (ppm)'])
    adj_pval=len(assignedresults[assignedresults[filterby]==c])

    #Calculate statistic from (how similar distributions are)
    ks=ks_2samp(assignedresults[assignedresults[filterby]==c]['m/z Error (ppm)'],
                   assignedresults[assignedresults['Molecular class']==ref_class]['m/z Error (ppm)'])

    if adj_pval>4:
        if ranksum.pvalue*len(assignedresults['rank class'].unique())>0.05:
            if ks.statistic<0.3:
                assignedresults.loc[assignedresults[filterby]==c,'class flag']=1


results1=assignedresults[assignedresults['class flag']==1]
results2=assignedresults[assignedresults['class flag']==0]

print(len(results1))
print(len(results2))

#### Plot and save error distribution figure

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
fig.set_size_inches(8, 6)

sns.scatterplot(x='m/z',y='m/z Error (ppm)',hue='Molecular class',data=results1,ax=ax1, edgecolor='none')
ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
ax1.set_title('a', fontweight='bold', loc='left')

sns.kdeplot(x='m/z Error (ppm)',data=results1,hue='Molecular class',ax=ax2,legend=False,common_norm=False)
ax2.set_title('b', fontweight='bold', loc='left')

sns.scatterplot(x='m/z',y='m/z Error (ppm)',hue='Molecular class',data=results2,ax=ax3, edgecolor='none')
ax3.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
ax3.set_title('c', fontweight='bold', loc='left')

sns.kdeplot(x='m/z Error (ppm)',data=results2,hue='Molecular class',ax=ax4,legend=False)
ax4.set_title('d', fontweight='bold', loc='left')

fig.tight_layout()

fig.savefig(file_location+'CoreLCMS_FigS3_rev.eps',dpi=300,format='eps')
fig.savefig(file_location+'CoreLCMS_FigS3_rev.pdf',dpi=300,format='pdf')

plt.show()
# Calculate atomic stoichiometries and Nominal Oxidation State of Carbon (NOSC)
#results['O/C']=results['O']/results['C']
#results['H/C']=results['H']/results['C']
#results['N/C']=results['N']/results['C']
#results['NOSC'] =  4 -(4*results['C'] + results['H'] - 3*results['N'] - 2*results['O'])/results['C']

#print('All assignments:', len(results))