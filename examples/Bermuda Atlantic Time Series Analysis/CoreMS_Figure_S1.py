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
warnings.filterwarnings("ignore")
import sys
sys.path.append("./")

# Change the current working directory to where CoreMS is located
os.chdir('/Users/boiteaur/Desktop/CoreMS_metallomics/CoreMS/')

# Import required modules
from corems.mass_spectra.input import rawFileReader

###### Set file folder and THERMO RAW file name here:
file_location='/Users/boiteaur/Desktop/Major projects/Bermuda Atlantic Time Series data processing/Thermo RAW data'
file="RMB_190828_BATSpooled_30.RAW" #pooled sample for formula assignments
bfile="RMB_190828_BATS24_blnk.RAW" #blank sample

parser = rawFileReader.ImportMassSpectraThermoMSFileReader(file_location+'/'+file)
bparser = rawFileReader.ImportMassSpectraThermoMSFileReader(file_location+'/'+bfile)

tic=parser.get_tic(ms_type='MS')[0]
btic=bparser.get_tic(ms_type='MS')[0]

tic_df=pd.DataFrame({'Time': tic.time,'Intensity': tic.tic,'Sample':'Pooled NOM'})
btic_df=pd.DataFrame({'Time': btic.time,'Intensity': btic.tic,'Sample':'Blank'})

tic_df=pd.concat([tic_df,btic_df])
tic_df=tic_df.reset_index()

fig, (ax) = plt.subplots(1)
sns.lineplot(x='Time',y='Intensity',data=tic_df,ax=ax, hue='Sample')
ax.set_xlabel('Time (min)')
ax.set_ylabel('Total Ion Current Intensity')
ax.set_xlim(0,35)
ax.axvline(x=4,color='black')
ax.axvline(x=30,color='black')
#ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

fig.savefig(file_location+'CoreLCMS_FigS1.eps',dpi=300,format='eps')
fig.savefig(file_location+'CoreLCMS_FigS1.pdf',dpi=300,format='pdf')

plt.show()