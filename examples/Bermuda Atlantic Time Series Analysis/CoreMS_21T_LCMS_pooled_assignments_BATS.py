# Python script for testing molecular formula assignment criteria on a single pooled sample. 
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

# Import required modules
import matplotlib.pyplot as plt
from corems.mass_spectra.input import rawFileReader
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.encapsulation.factory.parameters import MSParameters
from corems.mass_spectrum.calc.Calibration import MzDomainCalibration

######## Set files here 
# Set file folder and THERMO RAW file name here:
file_location='/Users/boiteaur/Desktop/Major projects/Bermuda Atlantic Time Series data processing/Thermo RAW data/'
file="RMB_190828_BATSpooled_30.RAW" #pooled sample for formula assignments
refmasslist = file_location+"cal_pos.ref"

# Change the current working directory to where CoreMS is located
#os.chdir('/Users/boiteaur/Desktop/CoreMS_metallomics/CoreMS/')

### Set time bins in minutes
interval=2
timerange=[0,30]
internal_cal_setting='Y' # Should be 'Y' to perform internal calibration.

#Molecular search parameters. 
MSParameters.molecular_search.error_method = 'None'
MSParameters.molecular_search.min_ppm_error = -0.25
MSParameters.molecular_search.max_ppm_error = 0.25
MSParameters.molecular_search.ion_charge = 1

MSParameters.mass_spectrum.threshold_method = 'signal_noise'
MSParameters.mass_spectrum.s2n_threshold=3
MSParameters.ms_peak.peak_min_prominence_percent = 0.1

MSParameters.mass_spectrum.min_picking_mz=200
MSParameters.mass_spectrum.max_picking_mz=800

####### End of parameters

MSParameters.molecular_search.url_database = "postgresql+psycopg2://coremsappdb:coremsapppnnl@localhost:5432/coremsapp"
MSParameters.molecular_search.score_method = "prob_score"
MSParameters.molecular_search.output_score_method = "prob_score"

print("Loading file: "+ file)
# Read in sample list and load MS data
MSfiles={}
parser = rawFileReader.ImportMassSpectraThermoMSFileReader(file_location+file)

MSfiles[file]=parser

#Function to calibrate the spectra in an LCMS run
def lcmsspectra_cal(parser,interval,timerange,internal_cal_setting):
    
    tic=parser.get_tic(ms_type='MS')[0]
    tic_df=pd.DataFrame({'time': tic.time,'scan': tic.scans})

    times=list(range(timerange[0],timerange[1],interval))

    calibrated_spectra={}
    
    for timestart in times:
        print(timestart)
        #Retrieve TIC for MS1 scans over the time range between 'timestart' and 'timestop' 
        scans=tic_df[tic_df.time.between(timestart,timestart+interval)].scan.tolist()
        
        #Now, get an average mass spectrum and list the centroided m/z values of the spectrum. One of these should be the molecule of interest.
        mass_spectrum = parser.get_average_mass_spectrum_by_scanlist(scans)
            
        mass_spectrum.settings.calib_sn_threshold = 5
        mass_spectrum.settings.min_calib_ppm_error = -3
        mass_spectrum.settings.max_calib_ppm_error = 3
        
        if(internal_cal_setting=='Y'):
            MzDomainCalibration(mass_spectrum, refmasslist,mzsegment=[0,1000]).run()

        calibrated_spectra[timestart]=mass_spectrum
    
    return(calibrated_spectra)

#Function to build formula assignment lists from calibrated spectra
def lcmsformula(spectra_dict):
    for key in spectra_dict:    
        print(key)
        SearchMolecularFormulas(spectra_dict[key], first_hit=True).run_worker_mass_spectrum()
        spectra_dict[key].percentile_assigned(report_error=True)

MSspectra={}
for file in MSfiles:
    MSspectra[file]=lcmsspectra_cal(MSfiles[file],interval,timerange,internal_cal_setting)



# Core Molecular formula search

MSParameters.molecular_search.min_dbe = 0
MSParameters.molecular_search.max_dbe = 20

#First, evaluate wide attribution criteria
MSParameters.molecular_search.usedAtoms['C'] = (4,50)
MSParameters.molecular_search.usedAtoms['H'] = (4,100)
MSParameters.molecular_search.usedAtoms['O'] = (1,20)
MSParameters.molecular_search.usedAtoms['N'] = (0,8)
MSParameters.molecular_search.usedAtoms['S'] = (0,0)
MSParameters.molecular_search.usedAtoms['Cu'] = (0,0)
MSParameters.molecular_search.usedAtoms['Si'] = (0,0)
MSParameters.molecular_search.adduct_atoms_pos: tuple = ('Na')
MSParameters.molecular_search.isProtonated = True
MSParameters.molecular_search.isRadical = False
MSParameters.molecular_search.isAdduct = True
MSParameters.molecular_search.max_oc_filter=1.2
MSParameters.molecular_search.max_hc_filter=3
#MSParameters.molecular_search.adduct_atoms_pos: tuple = ('Na', 'K')
lcmsformula(MSspectra[file])

#First, evaluate wide attribution criteria
MSParameters.molecular_search.usedAtoms['C'] = (4,50)
MSParameters.molecular_search.usedAtoms['H'] = (4,100)
MSParameters.molecular_search.usedAtoms['O'] = (1,20)
MSParameters.molecular_search.usedAtoms['N'] = (0,0)
MSParameters.molecular_search.usedAtoms['S'] = (0,3)
MSParameters.molecular_search.usedAtoms['Cu'] = (0,0)
MSParameters.molecular_search.usedAtoms['Si'] = (0,0)
MSParameters.molecular_search.isProtonated = True
MSParameters.molecular_search.isRadical = False
MSParameters.molecular_search.isAdduct = True
MSParameters.molecular_search.max_oc_filter=2
MSParameters.molecular_search.max_hc_filter=6

lcmsformula(MSspectra[file])

#First, evaluate wide attribution criteria
MSParameters.molecular_search.usedAtoms['C'] = (4,50)
MSParameters.molecular_search.usedAtoms['H'] = (4,100)
MSParameters.molecular_search.usedAtoms['O'] = (1,20)
MSParameters.molecular_search.usedAtoms['N'] = (0,0)
MSParameters.molecular_search.usedAtoms['S'] = (0,0)
MSParameters.molecular_search.usedAtoms['Cu'] = (0,0)
MSParameters.molecular_search.usedAtoms['Si'] = (0,10)
MSParameters.molecular_search.isProtonated = True
MSParameters.molecular_search.isRadical = False
MSParameters.molecular_search.isAdduct = False
MSParameters.molecular_search.max_oc_filter=2
MSParameters.molecular_search.max_hc_filter=6

lcmsformula(MSspectra[file])

masterresults={}
for file in MSspectra:
    results=[]
    for timebin in MSspectra[file]:
        print(file)
        print(timebin)
        assignments=MSspectra[file][timebin].to_dataframe()
        assignments['Time']=timebin
        results.append(assignments)

    results=pd.concat(results,ignore_index=True)
    results['File']=file
    masterresults[file]=results

allresults=pd.concat(masterresults.values())


#Fill zeros for elements involved in ratio calcs. 
elements=['C','H','O','N']
for element in elements:
    allresults[element]=allresults[element].fillna(0)

allresults['Molecular class']=allresults['Molecular Formula'].str.replace('\d+', '').str.replace(' ', '')
allresults['Molecular class'][allresults['Heteroatom Class']=='unassigned']='unassigned'
allresults['Molecular class'][allresults['Is Isotopologue']==1]='Isotope'

results=allresults[allresults['Molecular class']!='unassigned']

# Calculate atomic stoichiometries and Nominal Oxidation State of Carbon (NOSC)
results['O/C']=results['O']/results['C']
results['H/C']=results['H']/results['C']
results['N/C']=results['N']/results['C']
results['NOSC'] =  4 -(4*results['C'] + results['H'] - 3*results['N'] - 2*results['O'])/results['C']

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

#fig.savefig(file_location+'Phycosphere_library_errorplot.eps',dpi=300,format='eps')
#fig.savefig(file_location+'Phycosphere_library_errorplot.pdf',dpi=300,format='pdf')


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

#Calculate Dispersity Index. 
EIC={}
for file in allresults['File'].unique():
    masses=allresults[allresults['File']==file]['m/z'].unique().tolist()
    EIC[file]=MSfiles[file].get_eics(target_mzs=masses,tic_data={},peak_detection=False,smooth=False)
    
dispersity=[]
for ind in allresults.index:
    current=allresults.loc[ind]
    time=[0,2]+current.Time
    file=current.File
    mass=current['m/z']
    chroma=pd.DataFrame({'EIC':EIC[file][0][mass].eic,'time':EIC[file][0][mass].time})
    chroma=chroma[chroma['time'].between(time[0],time[1])]
    chroma=chroma.sort_values(by='EIC',ascending=False)
    d=chroma[chroma.cumsum()['EIC']<0.5*chroma.sum()['EIC']].time.std()
    dispersity.append(d)

allresults['Dispersity']=dispersity

plt.show()

allresults.to_csv(file_location+'/'+'assigned_results.csv')