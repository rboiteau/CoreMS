# Python script for generating reference mass calibration file in CoreMS. 
# RMB  5/26/2023
# Contributors: Yuri Corilo, Will Kew, Christian Dewey, Rene Boiteau

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
#from corems.molecular_id.factory.classification import HeteroatomsClassification, Labels
#from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.encapsulation.factory.parameters import MSParameters
from corems.mass_spectrum.calc.Calibration import MzDomainCalibration
#from corems.molecular_id.factory.molecularSQL import MolForm_SQL

######## Set files here 
# Set file folder and THERMO RAW file name here:
file_location='/Users/boiteaur/Desktop/Lignin/'
file="AS_230524_lignin_only_1.raw" #pooled sample for formula assignments
refmasslist = file_location+"cal_pos.ref"

# Change the current working directory to where CoreMS is located
os.chdir('/Users/boiteaur/Desktop/CoreMS_metallomics/CoreMS/')

### Set time bins in minutes
interval=4
timerange=[0,28]
internal_cal_setting='Y' # Should be 'Y' to do internal calibration first, 'N' to skip internal calibration.


#Molecular search parameters. 
MSParameters.molecular_search.error_method = 'None'
MSParameters.molecular_search.min_ppm_error = -1
MSParameters.molecular_search.max_ppm_error = 1
MSParameters.molecular_search.ion_charge = 1

MSParameters.mass_spectrum.threshold_method = 'minima'
MSParameters.mass_spectrum.noise_threshold_std = 50
#MSParameters.mass_spectrum.threshold_method = 's2n'
#MSParameters.mass_spectrum.s2n_threshold=20
MSParameters.ms_peak.peak_min_prominence_percent = 0.1

MSParameters.mass_spectrum.min_picking_mz=100
MSParameters.mass_spectrum.max_picking_mz=650


# Core Molecular formula search
MSParameters.molecular_search.min_dbe = 0
MSParameters.molecular_search.max_dbe = 20

MSParameters.molecular_search.usedAtoms['C'] = (4,40)
MSParameters.molecular_search.usedAtoms['H'] = (4,80)
MSParameters.molecular_search.usedAtoms['O'] = (1,18)
MSParameters.molecular_search.usedAtoms['N'] = (0,1)
MSParameters.molecular_search.usedAtoms['S'] = (0,0)
MSParameters.molecular_search.usedAtoms['Si'] = (0,0)
MSParameters.molecular_search.isProtonated = True
MSParameters.molecular_search.isAdduct = False
MSParameters.molecular_search.max_oc_filter=1.2
MSParameters.molecular_search.max_hc_filter=3

####### End of parameters

MSParameters.molecular_search.url_database = "postgresql+psycopg2://coremsappdb:coremsapppnnl@localhost:5432/coremsapp"
MSParameters.molecular_search.score_method = "prob_score"
MSParameters.molecular_search.output_score_method = "prob_score"

print("Loading file: "+ file)
# Read in sample list and load MS data
MSfiles={}
parser = rawFileReader.ImportMassSpectraThermoMSFileReader(file_location+'/'+file)

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



#Here, we create a new reference mass list.
cal_list=results[results['Confidence Score']>.5]
cal_list=results[results['Ion Charge']==1]
cal_list=cal_list[cal_list['Molecular class']!='Isotope'].drop_duplicates(subset=['Molecular Formula'])

fig, (ax1) = plt.subplots(1,1)
sns.scatterplot(x='m/z',y='m/z Error (ppm)',hue='Molecular class',data=cal_list, ax=ax1, edgecolor='none')

plt.show()

cal=pd.DataFrame({'# Name':cal_list['Molecular Formula'], 'm/z value':cal_list['Calculated m/z'], 'charge':cal_list['Ion Charge'],' ion formula':cal_list['Molecular Formula'],'collision cross section [A^2]':cal_list['Ion Charge']})

cal.to_csv(refmasslist,sep='\t',index=False)