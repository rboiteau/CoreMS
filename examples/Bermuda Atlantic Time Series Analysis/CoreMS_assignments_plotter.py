# Import the os module
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

file_location='/Users/boiteaur/Desktop/Major projects/Bermuda Atlantic Time Series data processing/Thermo RAW data/'

allresults=pd.read_csv(file_location+'BATS_allfiles_assigned_results.csv')

allresults['Molecular class'][allresults['N']>4]=='CHON>4'
allresults['Molecular class'][allresults['S']>1]=='CHOS>1'

results=allresults[allresults['Molecular class']!='unassigned']

# Calculate atomic stoichiometries and Nominal Oxidation State of Carbon (NOSC)
results['O/C']=results['O']/results['C']
results['H/C']=results['H']/results['C']
results['N/C']=results['N']/results['C']
results['NOSC'] =  4 -(4*results['C'] + results['H'] - 3*results['N'] - 2*results['O'])/results['C']

results=allresults[allresults['Is Isotopologue']==0]
results=allresults[allresults['Molecular class']!='Unassigned']

results_1=results[(results['Molecular class'].isin(['CHO','CHON']))]
results_2=results[(~results['Molecular class'].isin(['CHO','CHON']))]

#results_2=results[(results['Molecular class']!='CHO') | (results['Molecular class']!='CHON<4')]

print('All assignments:', len(results))
print('Good assignments:', len(results_1))
print('Questionable assignments:', len(results_2))




fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
fig.set_size_inches(8, 5)


sns.scatterplot(x='m/z',y='m/z Error (ppm)',hue='Molecular class',data=results_1,ax=ax1, edgecolor='none')
ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
ax1.set_xlim(200,800)
sns.kdeplot(x='m/z Error (ppm)',data=results_1,hue='Molecular class',ax=ax2,legend=False)
ax2.set_xlim(-0.3,0.3)
fig.tight_layout()


sns.scatterplot(x='m/z',y='m/z Error (ppm)',hue='Molecular class',data=results_2,ax=ax3, edgecolor='none')
ax3.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
ax3.set_xlim(200,800)
sns.kdeplot(x='m/z Error (ppm)',data=results_2,hue='Molecular class',ax=ax4,legend=False)
ax4.set_xlim(-0.3,0.3)
fig.tight_layout()

fig.savefig(file_location+'CoreLCMS_FigS3.eps',dpi=300,format='eps')
fig.savefig(file_location+'CoreLCMS_FigS3.pdf',dpi=300,format='pdf')

assign_summary=[]

for time in allresults['Time'].unique():
    current={}
    current['Time']=time
    for mol_class in sorted(allresults['Molecular class'].unique()):
        current[mol_class]=len(allresults[(allresults['Molecular class']==mol_class) & (allresults['Time']==time)])
    assign_summary.append(current)
    #mzdiff=result['m/z'].sort_values(ascending=True).diff().iloc[1:]/result['m/z'].sort_values(ascending=True).iloc[1:]*1E6


df=pd.DataFrame(assign_summary)
df.plot.bar(x='Time',y=df.columns[1:],stacked=True,ylabel='Peaks')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)



plt.show()
print('All peaks:', len(allresults))

print('All monoisotopic assignments:', len(results))


#### Plot and save error distribution figure
fig, ((ax1, ax2)) = plt.subplots(1,2)
fig.set_size_inches(12, 6)

sns.scatterplot(x='m/z',y='m/z Error (ppm)',hue='Molecular class',data=results,ax=ax1, edgecolor='none')
ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
ax1.set_title('a', fontweight='bold', loc='left')

sns.kdeplot(x='m/z Error (ppm)',data=results,hue='Molecular class',ax=ax2,legend=False)
ax2.set_title('b', fontweight='bold', loc='left')

fig.tight_layout()

#fig.savefig(file_location+'library_errorplot.eps',dpi=300,format='eps')
#fig.savefig(file_location+'library_errorplot.pdf',dpi=300,format='pdf')


#### Plot library assignments over time

assign_summary=[]
for time in allresults['Time'].unique():
    current={}
    current['Time']=time
    for mol_class in allresults['Molecular class'].unique():
        current[mol_class]=len(allresults[(allresults['Molecular class']==mol_class) & (allresults['Time']==time)])
    assign_summary.append(current)
    #mzdiff=result['m/z'].sort_values(ascending=True).diff().iloc[1:]/result['m/z'].sort_values(ascending=True).iloc[1:]*1E6


df=pd.DataFrame(assign_summary)
df.plot.bar(x='Time',y=df.columns[1:],stacked=True,ylabel='Peaks')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)

plt.show()