# Python script for generating Figure S5: Extracted ion chromatogram (EIC) of cyanocobalamin internal standard across all samples. 
# RMB update 6/02/2023
# Contributors: Yuri Corilo, Will Kew, Christian Dewey, Rene Boiteau

##########
# Import the os module
import os
import pandas as pd
import numpy as np
import seaborn as sns
import warnings
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
import sys
sys.path.append("./")

# Change the current working directory to where CoreMS is located
os.chdir('/CoreMS')


### Set data directory and file names here

data_dir = '/CoreMS/usrdata/'
sample_list_name='BATS_sample_list.csv' #Sample list must contain column with header 'File'
featurelist_file='BATS_featurelist.csv'
clustered_featurelist_file='TableS1_featurelist.csv'

# Import required modules
from corems.mass_spectra.input import rawFileReader

#Read in sample list and load MS data
samplelist=pd.read_csv(data_dir+sample_list_name)

MSfiles={}
for file in samplelist['File'][samplelist['Sample type']=='sample']:
    parser = rawFileReader.ImportMassSpectraThermoMSFileReader(data_dir+file)
    MSfiles[file]=parser

samplelist=samplelist[samplelist['File'].isin(MSfiles.keys())]

### (1) Internal standard QC Check
### Module compares internal standard peak area across samples. Outliers are flagged.

stdmass=678.2918 # m/z of cycanocobalamin
std_timerange=[13.2,14] # retention time range of peak (min)

# Plot standard EIC over time range and integrate peak area. 
area=[]
rt=[]

fig, axs = plt.subplot_mosaic([['a','b']], figsize=(11,5), constrained_layout=True)
axs['a'].set(xlabel='Time (min)',ylabel='Intensity',title='Internal Standard EIC = '+str(stdmass) + ' m/z')

for file in MSfiles.keys():
    EIC=MSfiles[file].get_eics(target_mzs=[stdmass],tic_data={},peak_detection=False,smooth=False)
    df=pd.DataFrame({'EIC':EIC[0][stdmass].eic,'time':EIC[0][stdmass].time})
    df_sub=df[df['time'].between(std_timerange[0],std_timerange[1])]
    area.append(sum(df_sub['EIC']))
    rt.append(df_sub.time[df_sub.EIC==df_sub.EIC.max()].max())
    axs['a'].plot(df_sub['time'],df_sub['EIC']/1e7,label=file[11:])
axs['a'].legend(loc='center left', bbox_to_anchor=(1, 0.5))
axs['a'].set_title('a', fontweight='bold', loc='left')
axs['a'].set_ylabel('Intensity (x 1e7)')

samplelist['qc_area']=area
samplelist['qc_rt']=rt


# Flag outliers with peak area greater than 2x standard deviation of the mean 

peak_stdv=samplelist[samplelist['Sample type']=='sample'].qc_area.std()
peak_mean=samplelist[samplelist['Sample type']=='sample'].qc_area.mean()

samplelist['qc_pass']=0
for i in samplelist.index:
    if (abs(samplelist.qc_area[i]-peak_mean)<2*peak_stdv):
        samplelist.qc_pass[i]=1

print(str(samplelist[samplelist['Sample type']=='sample'].qc_pass.sum()) + ' pass of ' + str(len(samplelist[samplelist['Sample type']=='sample'])))

peak_stdv=samplelist[samplelist.qc_pass==1].qc_area.std()

print(str(round(peak_stdv/peak_mean*100,1))+' % std dev')

#Create plot of overlaid standard EICs
sns.histplot(x='qc_area',data=samplelist,ax=axs['b'])
axs['b'].set_xlabel('Internal Standard Peak Area')
axs['b'].set_xlim(0,20e7)
axs['b'].set_title('b', fontweight='bold', loc='left')

plt.savefig(data_dir+'CoreLCMS_FigS5.eps',dpi=300,format='eps')
plt.savefig(data_dir+'CoreLCMS_FigS5.pdf',dpi=300,format='pdf')

plt.show()