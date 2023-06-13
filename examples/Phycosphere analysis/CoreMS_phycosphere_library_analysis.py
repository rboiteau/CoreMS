# Python script for batch molecular formula assignments 
# RMB Last updated  5/30/2023
# Contributors: Yuri Corilo, Will Kew, Christian Dewey, Rene Boiteau

# Applied to Bermuda Atlantic Time Series pooled sample 
##########

# Import the os module
import os
import pandas as pd
import numpy as np
import seaborn as sns
import warnings
warnings.filterwarnings("ignore")
import sys
sys.path.append("./")
from pathlib import Path

# Change the current working directory to where CoreMS is located
#os.chdir('/Users/boiteaur/Desktop/CoreMS_metallomics/CoreMS/')

# Import required modules
import matplotlib.pyplot as plt
from corems.mass_spectra.input import rawFileReader
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.encapsulation.factory.parameters import MSParameters
from corems.mass_spectrum.calc.Calibration import MzDomainCalibration

######## Set files here 
# Set file folder and THERMO RAW file name here:
file_location='/Users/boiteaur/Desktop/Major projects/Phycosphere LCMS/'
sample_list_name='Phycosphere_samplelist.csv' #Sample list must contain column with header 'File'
refmasslist = file_location+"siloxanes_pos.ref"
allresults_file='Phycosphere_21T_pooled_formula_library.csv'
library_file='Phycosphere_21T_pooled_formula_unique_library.csv'
#featurelist_file='AZTDI_alignedfeaturelist1_gapfilled.csv'
featurelist_file='Phycosphere_Height_0_2023661719.csv'
m_rt='Average Rt(min)'
m_mz='Average Mz'
offset=-0.4 #rt of library - retention time of files 
threshold=0.01 #Mass accuracy of metabolomic data in Da. 

MSParameters.mass_spectrum.min_picking_mz=150
MSParameters.mass_spectrum.max_picking_mz=800
MSParameters.mass_spectrum.threshold_method = 'log'
MSParameters.mass_spectrum.log_nsigma=400
MSParameters.ms_peak.peak_min_prominence_percent = 0.1

#Load MS data from sample list as MSfiles dictionary (keys=file name, values= parser objects)
samplelist=pd.read_csv(file_location+sample_list_name)
MSfiles={}
#for file in samplelist['File']:
    #parser = rawFileReader.ImportMassSpectraThermoMSFileReader(file_location+file)
    #MSfiles[file]=parser

allresults=pd.read_csv(file_location+allresults_file)


#Fill zeros for elements involved in ratio calcs. 
elements=['C','H','O','N','P','S']

#Annotate results based on stoichiometries
allresults['Stoichiometric classification']=0

for element in elements:
    if element in list(allresults.columns):
        allresults[element]=allresults[element].fillna(0)
    else:
        allresults[element]=0

allresults=allresults[allresults['N']<5]


assignedresults=allresults[allresults['Is Isotopologue']==0]
#assignedresults=allresults[allresults['C']<1]


# Calculate atomic stoichiometries and Nominal Oxidation State of Carbon (NOSC)
assignedresults['O/C']=assignedresults['O']/assignedresults['C']
assignedresults['H/C']=assignedresults['H']/assignedresults['C']
assignedresults['N/C']=assignedresults['N']/assignedresults['C']
assignedresults['P/C']=assignedresults['P']/assignedresults['C']
assignedresults['N/P']=assignedresults['N']/assignedresults['P']
assignedresults.loc[assignedresults['P']==0,'N/P']=0

assignedresults['NOSC'] =  4 -(4*assignedresults['C'] + assignedresults['H'] - 3*assignedresults['N'] - 2*assignedresults['O'])/assignedresults['C']

assignedresults['Stoichiometric classification']='Unclassified'

assignedresults.loc[(assignedresults['O/C']<=0.6) & 
                    (assignedresults['H/C']>=1.32) & 
                    (assignedresults['N/C']<=0.126) &
                    (assignedresults['P/C']<0.35)
                    ,'Stoichiometric classification'] = 'Lipid'

assignedresults.loc[(assignedresults['O/C']<=0.6) & 
                    (assignedresults['H/C']>=1.32) & 
                    (assignedresults['N/C']<=0.126) &
                    (assignedresults['P/C']<0.35) &
                    (assignedresults['P']>0)
                    ,'Stoichiometric classification'] = 'Phospholipid'

assignedresults.loc[(assignedresults['O/C']>=0.61) & 
                    (assignedresults['H/C']>=1.45) & 
                    (assignedresults['N/C']>0.07) & 
                    (assignedresults['N/C']<=0.2) & 
                    (assignedresults['P/C']<0.3) & 
                    (assignedresults['O']>=3) &
                    (assignedresults['N']>=1)
                    ,'Stoichiometric classification'] = 'A-Sugars'

assignedresults.loc[(assignedresults['O/C']>=0.8) & 
                    (assignedresults['H/C']>=1.65) & 
                    (assignedresults['H/C']<2.7) &
                    (assignedresults['O']>=3) &
                    (assignedresults['N']==0)
                    ,'Stoichiometric classification'] = 'Carbohydrates'

assignedresults.loc[(assignedresults['O/C']>=0.5) & 
                    (assignedresults['O/C']<1.7) & 
                    (assignedresults['H/C']>1) & 
                    (assignedresults['H/C']<1.8) &
                    (assignedresults['N/C']>=0.2) & 
                    (assignedresults['N/C']<=0.5) & 
                    (assignedresults['N']>=2) &
                    (assignedresults['P']>=1) &
                    (assignedresults['S']==0) &
                    (assignedresults['Calculated m/z']>305) &
                    (assignedresults['Calculated m/z']<523)
                    ,'Stoichiometric classification'] = 'Nucleotides'

assignedresults.loc[(assignedresults['O/C']<=1.15) & 
                    (assignedresults['H/C']<1.32) & 
                    (assignedresults['N/C']<0.126) &
                    (assignedresults['P/C']<=0.2) 
                    ,'Stoichiometric classification'] = 'Phytochemicals'

assignedresults.loc[(assignedresults['S']>0)
                    ,'Stoichiometric classification'] = 'Organosulfur'

assignedresults.loc[(assignedresults['O/C']>0.12) & 
                    (assignedresults['O/C']<=0.6) & 
                    (assignedresults['H/C']>0.9) & 
                    (assignedresults['H/C']<2.5) & 
                    (assignedresults['N/C']>=0.126) & 
                    (assignedresults['N/C']<=0.7) & 
                    (assignedresults['P/C']<0.17) & 
                    (assignedresults['N']>=1)
                    ,'Stoichiometric classification'] = 'Protein'

assignedresults.loc[(assignedresults['O/C']>0.6) & 
                    (assignedresults['O/C']<=1) & 
                    (assignedresults['H/C']>1.2) & 
                    (assignedresults['H/C']<2.5) & 
                    (assignedresults['N/C']>=0.2) & 
                    (assignedresults['N/C']<=0.7) & 
                    (assignedresults['P/C']<0.17) & 
                    (assignedresults['N']>=1)
                    ,'Stoichiometric classification'] = 'Protein'

#Create a list of all unique features and their abundances across samples.
uniquelist=[]
for time in assignedresults.Time.unique():
    current=assignedresults[assignedresults.Time==time]
    current=current.sort_values(by=['Peak Height'],ascending=False)
    currentunique=current.drop_duplicates(subset=['Molecular Formula'])
    currentunique=currentunique.set_index(['Molecular Formula'],drop=False)
    for file in assignedresults['File'].unique():
        current_file=current[current['File']==file]
        current_file=current_file.rename(columns={'Peak Height':file})
        current_file=current_file.set_index(['Molecular Formula'],drop=False)
        currentunique=currentunique.join(current_file[file])
    uniquelist.append(currentunique)

uniqueresults=pd.concat(uniquelist,ignore_index=True)

uniqueresults['Dispersity']=1

'''
EIC={}
for file in uniqueresults['File'].unique():
    masses=uniqueresults[uniqueresults['File']==file]['m/z'].unique().tolist()
    EIC[file]=MSfiles[file].get_eics(target_mzs=masses,tic_data={},peak_detection=False,smooth=False)


#Calculate Dispersity Index. 
dispersity=[]
for ind in uniqueresults.index:
    current=uniqueresults.loc[ind]
    time=[0,2]+current.Time
    file=current.File
    mass=current['m/z']
    chroma=pd.DataFrame({'EIC':EIC[file][0][mass].eic,'time':EIC[file][0][mass].time})
    chroma=chroma[chroma['time'].between(time[0],time[1])]
    chroma=chroma.sort_values(by='EIC',ascending=False)
    d=chroma[chroma.cumsum()['EIC']<0.5*chroma.sum()['EIC']].time.std()
    dispersity.append(d)

uniqueresults['Dispersity']=dispersity
'''

print('All peaks:', len(allresults))

print('All monoisotopic assignments:', len(assignedresults))

print('All unique assignments:', len(uniqueresults))

uniqueresults.to_csv(file_location+library_file)

#### Plot and save error distribution figure
fig, ((ax1, ax2)) = plt.subplots(1,2)
fig.set_size_inches(12, 6)

sns.scatterplot(x='m/z',y='m/z Error (ppm)',hue='Molecular class',data=uniqueresults,ax=ax1, edgecolor='none')
ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
ax1.set_title('a', fontweight='bold', loc='left')

sns.kdeplot(x='m/z Error (ppm)',data=uniqueresults,hue='Molecular class',ax=ax2,legend=False)
ax2.set_title('b', fontweight='bold', loc='left')

fig.tight_layout()

fig.savefig(file_location+'Phycosphere_library_errorplot.eps',dpi=300,format='eps')
fig.savefig(file_location+'Phycosphere_library_errorplot.pdf',dpi=300,format='pdf')

### Use masterresults library to annotate MZmine3 feature list
featurelist=pd.read_csv(file_location+'/'+ featurelist_file)

timebins=uniqueresults.Time.unique()

feature_annotations=[]
for i in featurelist.iterrows():
    current=i[1].to_dict()
    ctime=current[m_rt]+offset
    cmass=current[m_mz]
    match=(timebins-ctime)
    match=round(match[match<1].max()+ctime)

    annotations=allresults[(allresults['Time']==match) & (abs(allresults['m/z']-cmass)<threshold)]
    current['all library hits']=len(annotations)
    annotations=uniqueresults[(uniqueresults['Time']==match) & (abs(uniqueresults['m/z']-cmass)<threshold)]
    current['annotated library hits']=len(annotations)

    if len(annotations)>0:
        if len(annotations)>1:
            annotations=annotations[annotations['Peak Height']==max(annotations['Peak Height'])]
        current['theor m/z']=annotations['Calculated m/z'].to_numpy()[0]
        current['Molecular Formula']=annotations['Molecular Formula'].to_numpy()[0]
        current['Library Time']=annotations['Time'].to_numpy()[0]
        current['Library m/z error']=annotations['m/z Error (ppm)'].to_numpy()[0]
        current['Molecular class']=annotations['Molecular class'].to_numpy()[0]
        current['Library intensity']=annotations['Peak Height'].to_numpy()[0]
        current['Library ion charge']=annotations['Ion Charge'].to_numpy()[0]
        current['Library is isotopologue']=annotations['Is Isotopologue'].to_numpy()[0]
        current['m/z error']=(annotations['Calculated m/z'].to_numpy()[0]-cmass)/cmass*1e6
        current['O/C']=annotations['O/C'].to_numpy()[0]
        current['H/C']=annotations['H/C'].to_numpy()[0]
        current['N/C']=annotations['N/C'].to_numpy()[0]
        current['P/C']=annotations['P/C'].to_numpy()[0]
        current['DBE']=annotations['DBE'].to_numpy()[0]
        current['NOSC']=annotations['NOSC'].to_numpy()[0]
        current['Dispersity']=annotations['Dispersity'].to_numpy()[0]
        current['Stoichiometric classification']=annotations['Stoichiometric classification'].to_numpy()[0]

        for element in elements:
            current[element]=annotations[element].to_numpy()

    feature_annotations.append(current)


featurelist_annotated=pd.DataFrame(feature_annotations)

featurelist_annotated.to_csv(file_location+'/'+featurelist_file.replace('.csv','_annotated.csv'))


print(len(featurelist_annotated))
print(len(featurelist_annotated[featurelist_annotated['all library hits']>0]))

print(len(featurelist_annotated[featurelist_annotated['theor m/z']>0]))

#### Plot and save error distribution figure
fig, ((ax1, ax2)) = plt.subplots(1,2)
fig.set_size_inches(12, 6)

sns.scatterplot(x=m_rt,y='O/C',hue='Stoichiometric classification',data=featurelist_annotated,ax=ax1, edgecolor='none')
#ax1.set_xlim(200,800)
ax1.set_title('a', fontweight='bold', loc='left')
ax1.set_xlabel('Retention time (min)')

sns.scatterplot(x=m_rt,y='N/C',data=featurelist_annotated,hue='Stoichiometric classification',ax=ax2,legend=False)
#ax2.set_xlim(-0.3,0.3)
ax2.set_title('b', fontweight='bold', loc='left')
ax2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
ax2.set_xlabel('Retention time (min)')

fig.tight_layout()


#### Plot and save error distribution figure
fig, ((ax1, ax2)) = plt.subplots(1,2)
fig.set_size_inches(12, 6)

sns.scatterplot(x=m_rt,y='theor m/z',hue='Stoichiometric classification',data=featurelist_annotated,ax=ax1, edgecolor='none')
ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
#ax1.set_xlim(200,800)
ax1.set_title('a', fontweight='bold', loc='left')
ax1.set_xlabel('Retention time (min)')

sns.scatterplot(x=m_rt,y='DBE',data=featurelist_annotated,hue='Stoichiometric classification',ax=ax2,legend=False)
#ax2.set_xlim(-0.3,0.3)
ax2.set_title('b', fontweight='bold', loc='left')
ax2.set_xlabel('Retention time (min)')

fig.tight_layout()

plt.show()

