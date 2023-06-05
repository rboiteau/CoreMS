# Python script for generating Figure S2: Resolving power (m/z / Δm/z) required to separate DOM molecular peaks.
# RMB update 6/02/2023
# Contributors: Yuri Corilo, Will Kew, Christian Dewey, Rene Boiteau

##########

# Import the os module
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

file_location='/Users/boiteaur/Desktop/Major projects/Bermuda Atlantic Time Series data processing/Thermo RAW data/'

allresults=pd.read_csv(file_location+'BATS_allfiles_assigned_results.csv')

#sample_list_name='BATS_sample_list.csv' #Sample list must contain column with header 'File'

#Load MS data from sample list as MSfiles dictionary (keys=file name, values= parser objects)
#samplelist=pd.read_csv(file_location+sample_list_name)
#MSfiles={}
#for file in samplelist['File']:
#    parser = rawFileReader.ImportMassSpectraThermoMSFileReader(file_location+file)
#    MSfiles[file]=parser

#Calculate required resolving power and plot over time
res_summary=[]
mzinterval=100
mzint=range(200,900,mzinterval)
for time in allresults['Time'].unique():
    
    result=allresults[allresults['Time']==time]

    for file in allresults['File'].unique():

        result_sub=result[result['File']==file]

        for mz in mzint:
            result_sub_mz=result_sub[(result_sub['m/z']>mz) & (result_sub['m/z'] < mz + mzinterval)]
            mzvalues=result_sub_mz['m/z'].sort_values(ascending=True)
            differences=mzvalues.diff()
            #Resolve from peaks on either side
            mzdiff=pd.DataFrame({'left':differences[1:-1].to_list(),'right':differences[2:].to_list()})
            mzdiff=mzdiff.min(axis=1)
            mzvalues=mzvalues.iloc[1:-1].reset_index()['m/z']
            #mzdiff_res=mzvalues/mzdiff*2
            mzdiff_res=mzvalues/mzdiff*2
            for i in mzdiff_res.index:
                res_summary.append({'resolution':mzdiff_res[i],'mass':mzvalues[i],'Mass Range':str(mz)+'-'+str(mz+mzinterval) + ' m/z','file':file,'time':time})

res_summary_df=pd.DataFrame(res_summary)
fig, (ax) = plt.subplots(1)
sns.ecdfplot(x='resolution',hue='Mass Range',data=res_summary_df,ax=ax)
ax.set_xlim(0,800000)
ax.set_xlabel('Resolving power (m/z / '+r'$\Delta$'+'m/z)')
ax.set_ylabel('Fraction of peaks resolved')

fig.savefig(file_location+'CoreLCMS_FigS2.eps',dpi=300,format='eps')
fig.savefig(file_location+'CoreLCMS_FigS2.pdf',dpi=300,format='pdf')

plt.show()
