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
file_location='/Users/boiteaur/Desktop/Major projects/Bermuda Atlantic Time Series data processing/Thermo RAW data/'
sample_list_name='BATS_sample_list.csv' #Sample list must contain column with header 'File'
refmasslist = file_location+"Seawater_NOM_pos_recal3.ref"
savefile='BATS_allfiles_assigned_results_round4.csv'
dfiletype='.raw'

### Set time bins in minutes
interval=2
timerange=[2,36]
internal_cal_setting='Y' # Should be 'Y' to perform internal calibration.

#Molecular search parameters. 
MSParameters.molecular_search.error_method = 'None'
MSParameters.molecular_search.min_ppm_error = -0.4
MSParameters.molecular_search.max_ppm_error = 0.4
MSParameters.molecular_search.ion_charge = 1

MSParameters.mass_spectrum.min_calib_ppm_error = 0.5
MSParameters.mass_spectrum.max_calib_ppm_error = 2
MSParameters.mass_spectrum.calib_pol_order = 2
MSParameters.mass_spectrum.calib_sn_threshold = 7
MSParameters.mass_spectrum.min_picking_mz=200
MSParameters.mass_spectrum.max_picking_mz=900
MSParameters.mass_spectrum.threshold_method = 'signal_noise'
MSParameters.mass_spectrum.s2n_threshold=3

MSParameters.ms_peak.peak_min_prominence_percent = 0.001



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
        mass_spectrum.molecular_search_settings.usedAtoms['O'] = (1, 20)
        mass_spectrum.molecular_search_settings.usedAtoms['N'] = (0, 8)
        mass_spectrum.molecular_search_settings.usedAtoms['S'] = (0, 2)
        mass_spectrum.molecular_search_settings.usedAtoms['P'] = (0, 1)
        mass_spectrum.molecular_search_settings.usedAtoms['Si'] = (0, 0)
        mass_spectrum.molecular_search_settings.usedAtoms['Na'] = (0, 1)
        mass_spectrum.molecular_search_settings.isProtonated = True
        mass_spectrum.molecular_search_settings.isRadical = False
        mass_spectrum.molecular_search_settings.isAdduct = False
        mass_spectrum.molecular_search_settings.max_oc_filter=1.2
        mass_spectrum.molecular_search_settings.max_hc_filter=3
    
        SearchMolecularFormulas(mass_spectrum, first_hit=True).run_worker_mass_spectrum()
        mass_spectrum.percentile_assigned(report_error=True)
        #Add assigned spectrum to the calibrated_spectra dictionary
        calibrated_spectra[timestart]=mass_spectrum

    return(calibrated_spectra)

masterresults={}

for file in MSfiles:
    print(file)
    MSspectra=lcms_cal_assign(MSfiles[file],interval,timerange,internal_cal_setting)
    results=[]
    for timebin in MSspectra:
        print(timebin)
        assignments=MSspectra[timebin].to_dataframe()
        assignments['Time']=timebin
        results.append(assignments)
    allresults=pd.concat(results,ignore_index=True)
    allresults['File']=file
    allresults['Molecular class']=allresults['Molecular Formula'].str.replace('\d+', '').str.replace(' ', '')
    allresults['Molecular class'][allresults['Heteroatom Class']=='unassigned']='unassigned'
    allresults['Molecular class'][allresults['Is Isotopologue']==1]='Isotope'

    masterresults[file]=allresults

    allresults.to_csv(file_location+file.replace(dfiletype,'_mf_assignments_round4.csv'))

    assignedresults=allresults[allresults['Is Isotopologue']==0]

    #### Plot and save error distribution figure
    fig, ((ax1, ax2)) = plt.subplots(1,2)
    fig.set_size_inches(10, 5)

    sns.scatterplot(x='m/z',y='m/z Error (ppm)',hue='Molecular class',data=assignedresults,ax=ax1, edgecolor='none')
    ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)
    ax1.set_title('a', fontweight='bold', loc='left')

    sns.kdeplot(x='m/z Error (ppm)',data=assignedresults,hue='Molecular class',ax=ax2,legend=False)
    ax2.set_title('b', fontweight='bold', loc='left')

    fig.tight_layout()

    fig.savefig(file_location+file.replace(dfiletype,'_errorplot_round4.pdf'),dpi=300,format='pdf')

allresults=pd.concat(masterresults.values())

print('All peaks:', len(allresults))

print('All monoisotopic assignments:', len(allresults))

allresults.to_csv(file_location+savefile)

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

fig.savefig(file_location+savefile.replace('.csv','_assigment_plot_round4.eps'),dpi=300,format='eps')

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

