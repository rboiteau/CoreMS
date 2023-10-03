# Python script for generating Figure 5: Depth profiles and Extracted ion chromatograms of representative features. 
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


clusteredlist=clusteredlist[clusteredlist['File'].isin(samplelist['File'])]

depth=[]
for file in clusteredlist['File'].unique():
    d=samplelist[samplelist['File']==file]['Depth'].iloc[0]
    depth.append(d)


# Import required modules
from corems.mass_spectra.input import rawFileReader

#Load MS data from sample list as MSfiles dictionary (keys=file name, values= parser objects)

files=['RMB_190828_BATS18_100m.raw', 'RMB_190828_BATS19_800m.raw', 'RMB_190828_BATS20_5m.raw']


MSfiles={}
for file in files:
    parser = rawFileReader.ImportMassSpectraThermoMSFileReader(data_dir+file)
    MSfiles[file]=parser



### Module plots examples of individual feature profiles and EICs

matplotlib.rcParams['pdf.fonttype'] = '42'

def profile_eic_plot(mf,time):
    result=clusteredlist[(clusteredlist['Molecular Formula']==mf) & (clusteredlist['Time']==time)].squeeze()
    intensity=result[clusteredlist['File'].unique()].fillna(0)/1e5
    fig, (ax1, ax2) = plt.subplots(1,2)
    fig.set_size_inches(8,3.5)
    ax1.scatter(x=intensity,y=depth)
    ax1.set(xlabel='Average Intensity (x1e5)')
    ax1.set(xlim=[0,max(intensity)*1.05])
    ax1.set(ylabel='Depth (m)')
    ax1.set(title=result['Molecular Formula']+' Time:' + str(result['Time']) +' '+ result['cluster'].replace('\d+', ''))
    #ax1.axes.get_xaxis().set_major_locator(MaxNLocator(integer=True))
    ax1.invert_yaxis()

    stdmass=result['m/z']
    std_timerange=[2,36]

    area=[]
    rt=[]
    ax2.set(xlabel='Time (min)',ylabel='Intensity (x1e5)',title='EIC= '+str(stdmass.round(4))+' $\it{m/z}$')
    ax2.axes.get_yaxis().set_major_locator(ticker.MaxNLocator(integer=True))
    for file in files:

        EIC=MSfiles[file].get_eics(target_mzs=[stdmass],tic_data={},peak_detection=False,smooth=False)
        df=pd.DataFrame({'EIC':EIC[0][stdmass].eic,'time':EIC[0][stdmass].time})
        df_sub=df[df['time'].between(std_timerange[0],std_timerange[1])]
        area.append(sum(df_sub['EIC']))
        rt.append(df_sub.time[df_sub.EIC==df_sub.EIC.max()].max())
        ax2.plot(df_sub['time'],df_sub['EIC']/1e5,label=file[11:])

    ax2.legend(loc='upper left')
    fig.tight_layout()
    fig.savefig(data_dir+'CoreLCMS_Fig5'+mf+'.pdf',dpi=300,format='pdf')
    
profile_eic_plot('C30 H55 O7 N1',24)
profile_eic_plot('C29 H47 O7 N1',24)
profile_eic_plot('C29 H41 O7 N1',24)
profile_eic_plot('C15 H24 O6',24)

