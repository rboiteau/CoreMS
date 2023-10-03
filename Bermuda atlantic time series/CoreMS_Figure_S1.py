# Python script for generating Figure S1: Total ion chromatogram of the pooled sample and the C18 blank. 
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
warnings.filterwarnings('ignore')
import sys
sys.path.append('./')
os.chdir('/CoreMS')

# Import required modules
from corems.mass_spectra.input import rawFileReader

###### Set file folder and THERMO RAW file name here:
data_dir = '/CoreMS/usrdata/'
file="RMB_190828_BATSpooled_12.RAW" # Sample data file
file2="RMB_190828_BATSpooled_21.RAW" # Sample data file
file3="RMB_190828_BATSpooled_30.RAW" # Sample data file
bfile="RMB_190828_BATS24_blnk.RAW" # Blank data file

###### End user input


parser = rawFileReader.ImportMassSpectraThermoMSFileReader(data_dir+file)
tic=parser.get_tic(ms_type='MS')[0]
tic_df=pd.DataFrame({'Time': tic.time,'Intensity': tic.tic,'Sample':'Pooled Sample 1'})
df=tic_df


parser = rawFileReader.ImportMassSpectraThermoMSFileReader(data_dir+file2)
tic=parser.get_tic(ms_type='MS')[0]
tic_df=pd.DataFrame({'Time': tic.time,'Intensity': tic.tic,'Sample':'Pooled Sample 2'})
df=pd.concat([df,tic_df])
df=df.reset_index()

parser = rawFileReader.ImportMassSpectraThermoMSFileReader(data_dir+file3)
tic=parser.get_tic(ms_type='MS')[0]
tic_df=pd.DataFrame({'Time': tic.time,'Intensity': tic.tic,'Sample':'Pooled Sample 3'})
df=pd.concat([df,tic_df])
df=df.reset_index()

parser = rawFileReader.ImportMassSpectraThermoMSFileReader(data_dir+bfile)
tic=parser.get_tic(ms_type='MS')[0]
tic_df=pd.DataFrame({'Time': tic.time,'Intensity': tic.tic,'Sample':'Blank'})
df=pd.concat([df,tic_df])
#df=df.reset_index()



fig, (ax) = plt.subplots(1)
sns.lineplot(x='Time',y='Intensity',data=df,ax=ax, hue='Sample')
ax.set_xlabel('Time (min)')
ax.set_ylabel('Total Ion Current Intensity')
ax.set_xlim(0,36)
#ax.axvline(x=4,color='black')
#ax.axvline(x=30,color='black')
#ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

fig.savefig(data_dir+'CoreLCMS_FigS1.pdf',dpi=300,format='pdf')

plt.show()