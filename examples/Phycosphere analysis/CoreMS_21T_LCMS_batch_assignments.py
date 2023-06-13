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
savefile='Phycosphere_21T_pooled_formula_library.csv'
mstype='.RAW'
### Set time bins in minutes
interval=4
timerange=[0,32]
internal_cal_setting='Y' # Should be 'Y' to perform internal calibration.

#Molecular search parameters. 
MSParameters.molecular_search.error_method = 'None'
MSParameters.molecular_search.min_ppm_error = -0.25
MSParameters.molecular_search.max_ppm_error = 0.25
MSParameters.molecular_search.ion_charge = 1

MSParameters.mass_spectrum.min_calib_ppm_error = 0
MSParameters.mass_spectrum.max_calib_ppm_error = 1
MSParameters.mass_spectrum.calib_pol_order = 2
MSParameters.mass_spectrum.calib_sn_threshold = 10
MSParameters.mass_spectrum.min_picking_mz=150
MSParameters.mass_spectrum.max_picking_mz=800
MSParameters.mass_spectrum.threshold_method = 'log'
MSParameters.mass_spectrum.log_nsigma=400
#MSParameters.mass_spectrum.threshold_method = 'signal_noise'
#MSParameters.mass_spectrum.s2n_threshold=3
#MSParameters.mass_spectrum.threshold_method = 'minima'
#MSParameters.mass_spectrum.noise_threshold_std = 6

MSParameters.ms_peak.peak_min_prominence_percent = 0.1



####### End of settable parameters

MSParameters.molecular_search.url_database = "postgresql+psycopg2://coremsappdb:coremsapppnnl@localhost:5432/coremsapp"
MSParameters.molecular_search.score_method = "prob_score"
MSParameters.molecular_search.output_score_method = "prob_score"

#Load MS data from sample list as MSfiles dictionary (keys=file name, values= parser objects)
samplelist=pd.read_csv(file_location+sample_list_name)
MSfiles={}
for file in samplelist['File']:
    parser = rawFileReader.ImportMassSpectraThermoMSFileReader(file_location+file)
    MSfiles[file]=parser

#Function to calibrate and assign formula to the spectra in an LCMS run
def lcms_cal_assign(parser,interval,timerange,internal_cal_setting):
    
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

        #Calibrate spectrum based on reference mass list:
        if(internal_cal_setting=='Y'):

            MzDomainCalibration(mass_spectrum, refmasslist,mzsegment=[0,1000]).run()

        #Assign molecular formula based on specified elemental criteria

        #First assignment iteration (CHON with adducts)
        mass_spectrum.molecular_search_settings.min_dbe = 0
        mass_spectrum.molecular_search_settings.max_dbe = 20
        mass_spectrum.molecular_search_settings.adduct_atoms_pos = ('Na',) #Note, this must be a tuple
        mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1, 50)
        mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4, 100)
        mass_spectrum.molecular_search_settings.usedAtoms['O'] = (1, 15)
        mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0, 8)
        mass_spectrum.molecular_search_settings.usedAtoms['S'] = (0, 1)
        mass_spectrum.molecular_search_settings.usedAtoms['P'] = (0, 1)
        mass_spectrum.molecular_search_settings.usedAtoms['Si'] = (0, 0)
        mass_spectrum.molecular_search_settings.isProtonated = True
        mass_spectrum.molecular_search_settings.isRadical = False
        mass_spectrum.molecular_search_settings.isAdduct = True
        mass_spectrum.molecular_search_settings.max_oc_filter=1.2
        mass_spectrum.molecular_search_settings.max_hc_filter=3
    
        SearchMolecularFormulas(mass_spectrum, first_hit=True).run_worker_mass_spectrum()
        mass_spectrum.percentile_assigned(report_error=True)
        '''
        #Second assignment iteration (Sulfur containing)
        mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1, 50)
        mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4, 100)
        mass_spectrum.molecular_search_settings.usedAtoms['O'] = (1, 20)
        mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0, 0)
        mass_spectrum.molecular_search_settings.usedAtoms['S'] = (1, 3)
        mass_spectrum.molecular_search_settings.usedAtoms['P'] = (0, 0)
        mass_spectrum.molecular_search_settings.usedAtoms['Na'] = (0, 1)
        mass_spectrum.molecular_search_settings.usedAtoms['Si'] = (0, 0)
        mass_spectrum.molecular_search_settings.isProtonated = True
        mass_spectrum.molecular_search_settings.isRadical = False
        mass_spectrum.molecular_search_settings.isAdduct = False
        mass_spectrum.molecular_search_settings.max_oc_filter=1.2
        mass_spectrum.molecular_search_settings.max_hc_filter=3

        #SearchMolecularFormulas(mass_spectrum, first_hit=True).run_worker_mass_spectrum()
        #mass_spectrum.percentile_assigned(report_error=True)

        #Third assignment iteration (Organo-Phosphates)
        mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1, 50)
        mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4, 100)
        mass_spectrum.molecular_search_settings.usedAtoms['O'] = (1, 20)
        mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0, 8)
        mass_spectrum.molecular_search_settings.usedAtoms['S'] = (0, 0)
        mass_spectrum.molecular_search_settings.usedAtoms['P'] = (1, 1)
        mass_spectrum.molecular_search_settings.usedAtoms['Na'] = (0, 1)
        mass_spectrum.molecular_search_settings.usedAtoms['Si'] = (0, 0)
        mass_spectrum.molecular_search_settings.isProtonated = True
        mass_spectrum.molecular_search_settings.isRadical = False
        mass_spectrum.molecular_search_settings.isAdduct = False
        mass_spectrum.molecular_search_settings.max_oc_filter=1.2
        mass_spectrum.molecular_search_settings.max_hc_filter=3

        #SearchMolecularFormulas(mass_spectrum, first_hit=True).run_worker_mass_spectrum()
        #mass_spectrum.percentile_assigned(report_error=True)

        #Fourth assignment iteration (Siloxanes)
        mass_spectrum.molecular_search_settings.usedAtoms['C'] = (1, 40)
        mass_spectrum.molecular_search_settings.usedAtoms['H'] = (4, 80)
        mass_spectrum.molecular_search_settings.usedAtoms['O'] = (1, 20)
        mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0, 0)
        mass_spectrum.molecular_search_settings.usedAtoms['S'] = (0, 0)
        mass_spectrum.molecular_search_settings.usedAtoms['P'] = (0, 0)
        mass_spectrum.molecular_search_settings.usedAtoms['Si'] = (2, 10)
        mass_spectrum.molecular_search_settings.isProtonated = True
        mass_spectrum.molecular_search_settings.isRadical = False
        mass_spectrum.molecular_search_settings.isAdduct = False
        mass_spectrum.molecular_search_settings.max_oc_filter=2
        mass_spectrum.molecular_search_settings.max_hc_filter=6

        #SearchMolecularFormulas(mass_spectrum, first_hit=True).run_worker_mass_spectrum()
        #mass_spectrum.percentile_assigned(report_error=True)
        '''
        #Add assigned spectrum to the calibrated_spectra dictionary
        calibrated_spectra[timestart]=mass_spectrum

    return(calibrated_spectra)

MSspectra={}
for file in MSfiles:
    print(file)
    MSspectra[file]=lcms_cal_assign(MSfiles[file],interval,timerange,internal_cal_setting)

masterresults={}
for file in MSspectra:
    results=[]
    for timebin in MSspectra[file]:
        print(file)
        print(timebin)
        assignments=MSspectra[file][timebin].to_dataframe()
        assignments['Time']=timebin
        #pd.DataFrame(assignments).to_csv(file_location+file.replace(mstype,'_CoreMS_annotation.csv'))
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

assignedresults=allresults[allresults['Is Isotopologue']==0]

allresults.to_csv(file_location+savefile)

# Calculate atomic stoichiometries and Nominal Oxidation State of Carbon (NOSC)
assignedresults['O/C']=assignedresults['O']/assignedresults['C']
assignedresults['H/C']=assignedresults['H']/assignedresults['C']
assignedresults['N/C']=assignedresults['N']/assignedresults['C']
assignedresults['NOSC'] =  4 -(4*assignedresults['C'] + assignedresults['H'] - 3*assignedresults['N'] - 2*assignedresults['O'])/assignedresults['C']

print('All peaks:', len(allresults))

print('All monoisotopic assignments:', len(assignedresults))



#### Plot and save error distribution figure
fig, ((ax1, ax2)) = plt.subplots(1,2)
fig.set_size_inches(12, 6)

sns.scatterplot(x='m/z',y='m/z Error (ppm)',hue='File',data=assignedresults,ax=ax1, edgecolor='none')
ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
ax1.set_title('a', fontweight='bold', loc='left')

sns.scatterplot(x='O/C',y='H/C',data=assignedresults,hue='File',ax=ax2,legend=False)
ax2.set_title('b', fontweight='bold', loc='left')

fig.tight_layout()

fig.savefig(file_location+'Phycosphere_library_errorplot.eps',dpi=300,format='eps')
fig.savefig(file_location+'Phycosphere_library_errorplot.pdf',dpi=300,format='pdf')

#### Plot library assignments over time
assign_summary=[]
for time in allresults['Time'].unique():
    current={}
    current['Time']=time
    for mol_class in allresults['Molecular class'].unique():
        current[mol_class]=len(allresults[(allresults['Molecular class']==mol_class) & (allresults['Time']==time)])
    assign_summary.append(current)

df=pd.DataFrame(assign_summary)
df.plot.bar(x='Time',y=df.columns[1:],stacked=True,ylabel='Peaks')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=1.,frameon=True)

#Calculate Dispersity Index. 
'''
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
'''
plt.show()

