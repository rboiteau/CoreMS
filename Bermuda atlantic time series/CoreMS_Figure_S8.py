# Python script for generating Figure S8: Dispersity illustration
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
import matplotlib.ticker as ticker
import matplotlib.backends.backend_pdf

warnings.filterwarnings('ignore')
import sys
sys.path.append('./')
os.chdir('/CoreMS')


### Set data directory and file names here

data_dir = '/CoreMS/usrdata/'
sample_list_name='BATS_sample_list.csv' #Sample list must contain column with header 'File'
featurelist_file='BATS_featurelist.csv'
clustered_featurelist_file='TableS1_featurelist.csv'

### End user input


clusteredlist=pd.read_csv(data_dir+clustered_featurelist_file).fillna(0)
samplelist=pd.read_csv(data_dir+sample_list_name)
samplelist=samplelist[samplelist['Sample type']=='sample']
samplelist['Depth']=samplelist['Depth (m)']
samplelist=samplelist[samplelist['Depth']>0]


# Import required modules
from corems.mass_spectra.input import rawFileReader


#Calculate Dispersity Index. 
#EIC={}
#for file in clusteredlist['File'].unique():
#    masses=clusteredlist[clusteredlist['File']==file]['m/z'].unique().tolist()
#    EIC[file]=MSfiles[file].get_eics(target_mzs=masses,tic_data={},peak_detection=False,smooth=False)
    
matplotlib.rcParams['pdf.fonttype'] = '42'

#Define function that plots EIC of a feature and subset of peaks used for dispersity index calculation.
def dispersity_plotter(mf,time,clusteredlist):
    fig, (ax1, ax2) = plt.subplots(2,1,sharex=True,gridspec_kw=dict(height_ratios=[0.5,3]))
    fig.set_size_inches(6,3)    
    current=clusteredlist[(clusteredlist['Molecular Formula']==mf) & (clusteredlist['Time']==time)].squeeze()
    time=[0,2]+current.Time
    file=current.File
    mass=current['m/z']
    parser = rawFileReader.ImportMassSpectraThermoMSFileReader(data_dir+file)
    EIC=parser.get_eics(target_mzs=[mass],tic_data={},peak_detection=False,smooth=False)
    chroma=pd.DataFrame({'EIC':EIC[0][mass].eic,'time':EIC[0][mass].time})
    chrom_all=chroma[chroma['time'].between(0,8)]
    chroma=chroma[chroma['time'].between(time[0],time[1])]
    chroma=chroma.sort_values(by='EIC',ascending=False)
    chroma2=chroma[chroma.cumsum()['EIC']<0.5*chroma.sum()['EIC']]
    d=chroma[chroma.cumsum()['EIC']<0.5*chroma.sum()['EIC']].time.std()
    sns.boxplot(x='time',data=chroma2,ax=ax1,color='#E19E27')
    ax1.set_xlabel(None)
    ax1.set_title(current['Molecular Formula']+', Time:' + str(current['Time']) + ' min, Dispersity Index:' + str(current['Dispersity'].round(2)))
    sns.lineplot(x='time',y='EIC',data=chrom_all,ax=ax2,color='gray')
    sns.scatterplot(x='time',y='EIC',data=chroma,ax=ax2)
    sns.scatterplot(x='time',y='EIC',data=chroma2,ax=ax2)
    ax2.set_xlabel('Time (min)')
    #fig.savefig(data_dir+'CoreLCMS_FigS8_'+current['Molecular Formula']+'.eps',dpi=300,format='eps')
    fig.savefig(data_dir+'CoreLCMS_FigS8_'+current['Molecular Formula']+'.pdf',dpi=300,format='pdf')

#Plot examples of dispersity index calculations (Fig S8) for selected illustrative features
dispersity_plotter('C15 H27 O4 N1',4,clusteredlist)
dispersity_plotter('C18 H22 O7',4,clusteredlist)
