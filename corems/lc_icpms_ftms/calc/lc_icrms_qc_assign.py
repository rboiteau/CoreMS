
####Christian 09_30_2022 

import os
from matplotlib import style
import pandas as pd
import numpy as np

#Generate clusters of ms features across depth.
from sklearn.cluster import AgglomerativeClustering
import seaborn as sns
from matplotlib.colors import LogNorm, Normalize

import warnings

warnings.filterwarnings("ignore")
import sys
sys.path.append("./")
from pathlib import Path
import matplotlib.pyplot as plt

from corems.mass_spectra.input import rawFileReader
from corems.molecular_id.factory.classification import HeteroatomsClassification, Labels
from corems.molecular_id.search.priorityAssignment import OxygenPriorityAssignment
from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
from corems.encapsulation.factory.parameters import MSParameters
from corems.mass_spectrum.calc.Calibration import MzDomainCalibration


class lc_icr_assign:
    """
    Christian Dewey 
    Rene Boiteau

    October, 2022
    """

    def __init__(self, raw_files_location):

        self.raw_files_location = raw_files_location

        if raw_files_location[-1] != '/':
            
            self.savedir = raw_files_location + '/'

        else:

            self.savedir = raw_files_location            

        self._unique_run = False

        self._get_file_dict()


    def _get_file_dict(self):

        self._filelist=os.listdir(self.raw_files_location)

        self.master_data_holder={}

        self._raw_filelist = []

        for file in self._filelist:

            if '.raw' in file:

                parser = rawFileReader.ImportMassSpectraThermoMSFileReader(self.raw_files_location+'/'+file)

                self.master_data_holder[file]={'parser': parser}

                self._raw_filelist.append(file)


    def get_b12_eic(self, stdfile = None, plot_eic = True):

        MSParameters.mass_spectrum.threshold_method = 'signal_noise'
        MSParameters.mass_spectrum.s2n_threshold =6
        MSParameters.ms_peak.peak_min_prominence_percent = 0.1

        if stdfile == None:

            for f in self._filelist:

                if ('std' in f) or ('Std' in f):

                    parser = self.master_data_holder[f]['parser']

                    break

        else:

            parser = self.master_data_holder[stdfile]['parser']

        tic = parser.get_tic(ms_type = 'MS')[0]

        mass_spectrum = parser.get_average_mass_spectrum_by_scanlist(tic.scans)

        mz_exp_np = np.array(mass_spectrum.mz_exp)

        b12_mass = list(mz_exp_np[np.where((mz_exp_np > 678) & (mz_exp_np < 680))])
            
        stdEIC=parser.get_eics(target_mzs=b12_mass,tic_data={},peak_detection=False,smooth=False)

        _,ax = plt.subplots()

        max_std_intensity = 0

        dominant_mz = None
        
        for mz in b12_mass:
        
            ax.plot(stdEIC[0][mz].time, stdEIC[0][mz].eic, label='%.4f ' %mz)
        
            ax.legend(loc='best', frameon=False)

            if max_std_intensity < max(stdEIC[0][mz].eic):

                max_std_intensity = max(stdEIC[0][mz].eic)

                dominant_mz = mz

        print('max internal std intensity: %.2f' %max_std_intensity)
        print('dominant std mz: %.4f' %dominant_mz)

        self._stdmz = dominant_mz

        if plot_eic:
           
            plt.show()

        else:

            return ax
       

    def _qc_b12_plot(self,parser,std_timerange, stdmass, linefmt=None):

        ax=plt.subplot()
        ax.set(xlabel='Time (min)',ylabel='Intensity',title='Std EIC= %.4f' %(stdmass))

        EIC=parser.get_eics(target_mzs=[self._stdmz],tic_data={},peak_detection=False,smooth=False)
        df=pd.DataFrame({'EIC':EIC[0][self._stdmz].eic,'time':EIC[0][stdmass].time})
        df_sub=df[df['time'].between(std_timerange[0],std_timerange[1])]

        if linefmt:
            ax.plot(df_sub['time'],df_sub['EIC'],linefmt)
        else:
            ax.plot(df_sub['time'],df_sub['EIC'])
        
        area=sum(df_sub['EIC'])
        rt=df_sub.time[df_sub.EIC==df_sub.EIC.max()].max()

        return ax, area, rt


    def run_b12_qc(self, b12_peakrange, showplot = True, exlcude_std = True):

        plt.figure(figsize=(8,5))

        for file in self.master_data_holder.keys():

            parser = self.master_data_holder[file]['parser']

            if 'qH2O' in file:

                continue
            
            elif ('StdMix' in file) | ('stdmix' in file):
                
                ax, area, rt = self._qc_b12_plot(parser, b12_peakrange, self._stdmz, linefmt=':k')
            
            else:
                
                ax, area, rt = self._qc_b12_plot(parser, b12_peakrange, self._stdmz)

            self.master_data_holder[file]['b12_area'] = area
            self.master_data_holder[file]['b12_rt'] = rt
        
        ax.legend(labels = list(self.master_data_holder.keys()), loc='center left', prop={'size': 6}, bbox_to_anchor=(1, 0.5))

        fig = plt.gcf()
        fig.tight_layout()
        fig.savefig(self.savedir + 'b12_eics.pdf',bbox_inches = 'tight')

        if showplot:
            
            plt.show()

        if exlcude_std:
        
            areas = [self.master_data_holder[f]['b12_area'] for f in self._raw_filelist if ( ('qH2O' not in f) and ('StdMix' not in f) )]

        else:

            areas = [self.master_data_holder[f]['b12_area'] for f in self._raw_filelist if ( 'qH2O' not in f  )]
        
        rts = [self.master_data_holder[f]['b12_rt'] for f in self._raw_filelist if ('qH2O' not in f) ]

        area_stdev = np.std(areas)
        area_mean = np.mean(areas)

        rt_stdev = np.std(rts)
        rt_mean = np.mean(rts)
        
        fig, (ax, axrt) = plt.subplots(1, 2, figsize=(6,4))

        ax.boxplot(areas, labels=['positive mode'])
        ax.set_ylabel('b12 peak area')
        ax.text(x = .05  , y = 0.95 , s = '%.1f%% std dev' %(area_stdev/area_mean*100),transform=ax.transAxes, fontsize = 'x-small')
        axrt.boxplot(rts, labels=['positive mode'])
        axrt.set_ylabel('b12 r.t. (min)')
        axrt.text(x = .05  , y = 0.95 , s = '%.1f%% std dev' %(rt_stdev/rt_mean*100),transform=axrt.transAxes, fontsize = 'x-small')
        fig.tight_layout()
        fig.savefig(self.savedir +'b12_qcs.pdf')
        
        if showplot:
    
            plt.show()


    def _assign_formula(self, parser, interval, timerange, refmasslist=None,corder=2):
        #Function to build formula assignment lists
        #Retrieve TIC for MS1 scans over the time range between 'timestart' and 'timestop' 
        tic=parser.get_tic(ms_type='MS')[0]
        tic_df=pd.DataFrame({'time': tic.time,'scan': tic.scans})

        times=list(range(timerange[0],timerange[1],interval))

        results=[]
        
        for timestart in times:
            print('timestart: %s' %timestart )
            scans=tic_df[tic_df.time.between(timestart,timestart+interval)].scan.tolist()

            mass_spectrum = parser.get_average_mass_spectrum_by_scanlist(scans)    
            mass_spectrum.molecular_search_settings.ion_charge = 1
            #mass_spectrum.mass_spectrum.settings.calib_sn_threshold
            #mass_spectrum.mass_spectrum.settings.calib_pol_order
            #mass_spectrum.recalibrate_mass_spectrum(mass_spectrum, imzmeas, mzrefs, order=2)
            #MzDomainCalibration(mass_spectrum, ref_file_location).run()

            if refmasslist:
                mass_spectrum.settings.min_calib_ppm_error = 10
                mass_spectrum.settings.max_calib_ppm_error = -10
                calfn = MzDomainCalibration(mass_spectrum, refmasslist)
                ref_mass_list_fmt = calfn.load_ref_mass_list(refmasslist)

                imzmeas, mzrefs = calfn.find_calibration_points(mass_spectrum, ref_mass_list_fmt,
                                                            calib_ppm_error_threshold=(-1, 1),
                                                            calib_snr_threshold=3)

                calfn.recalibrate_mass_spectrum(mass_spectrum, imzmeas, mzrefs, order=corder)


            SearchMolecularFormulas(mass_spectrum, first_hit=False).run_worker_mass_spectrum()

            mass_spectrum.percentile_assigned(report_error=True)

            assignments=mass_spectrum.to_dataframe()

            assignments['Time']=timestart

            results.append(assignments)
        
        results=pd.concat(results,ignore_index=True)

        return(results)    


    def assign_formula(self, interval = None, timerange = None, refmasslist = None, calorder = 2):

        self.complete_results = {}
        ii = 1
        for file in self.master_data_holder:

            print('\n\n' + file)
            print("%s of %s files" %(ii, len(self.master_data_holder.keys())))

            results = self._assign_formula(self.master_data_holder[file]['parser'], interval, timerange, refmasslist=refmasslist, corder=calorder)
            results['file'] = file 
            self.master_data_holder[file]['results'] = results
            self.complete_results[file] = results

            ii = ii + 1 


    def plot_reqd_resolving_power(self, onlysamples=True):

        #Make plots showing required resolving power.
        if onlysamples:

            only_samples = {k: self.complete_results[k] for k in self.complete_results.keys() if ( 'stdmix' not in k) and ('blank' not in k)  }

            self.complete_results_df=pd.concat(only_samples.values())

        else:

            self.complete_results_df=pd.concat(self.complete_results.values())

        diff_summary=[]
        res_summary=[]

        for time in self.complete_results_df['Time'].unique():

            result=self.complete_results_df[self.complete_results_df['Time']==time]

            for file in self.complete_results_df['file'].unique():

                result_sub=result[result['file']==file]

                result_sub=result_sub['m/z'].sort_values(ascending=True)
                #mzdiff=result_sub['m/z'].sort_values(ascending=True).diff().iloc[1:]
                #mzdiff_res=result_sub['m/z'].iloc[1:]/mzdiff
                differences=result_sub.diff()

                #Resolve from peaks on either side:
                mzdiff=pd.DataFrame({'left':differences[1:-1].to_list(),'right':differences[2:].to_list(),'mz':result_sub[1:-1].to_list()})
                mzdiff_min=abs(mzdiff.min(axis=1))
                mzdiff_res=mzdiff['mz']/(mzdiff_min/2)

                diff_summary.extend(mzdiff_min.tolist())
                res_summary.extend(mzdiff_res[mzdiff_res<1e6].tolist())
                #res_summary.extend(mzdiff_res.tolist())

        res=list(range(10,int(max(res_summary)),1000))

        count=[]

        for i in res:

            count.append(len([element for element in res_summary if element<i])/len(res_summary))

        #print(len(diff_summary))
        #print(len([element for element in diff_summary if element<0.01]))

        fig, ax = plt.subplots()

        ax.plot(res,count)

        for i in [0.8,0.9,0.95,0.99]:

            current=abs(np.array(count)-i).tolist()

            print(res[current.index(min(current))])

            ax.axvline(res[current.index(min(current))],color='black')

        ax.set_xlabel('Resolving power (m/z / dm/z)')
        
        ax.set_ylabel('Fraction of peaks resolved')

        fig.tight_layout()

        fig.savefig(self.savedir +'b12_qcs.pdf')
        
        plt.show()


    def assess_all_results(self):

        only_samples = {k: self.complete_results[k] for k in self.complete_results.keys() if ( 'stdmix' not in k) }
        self.complete_results_df=pd.concat(only_samples.values())

        self.all_results=self.complete_results_df[(self.complete_results_df['m/z']<800) & (self.complete_results_df['S/N']>3)]

        self.all_results['N']=self.all_results['N'].fillna(0)
        self.all_results['O']=self.all_results['O'].fillna(0)
        #self.all_results['S']=self.all_results['S'].fillna(0)
        #self.all_results['P']=self.all_results['P'].fillna(0)
        #self.all_results['Na']=self.all_results['Na'].fillna(0)

        self.all_results['mol_class']='Unassigned'
        self.all_results['mol_class'][self.all_results['C']>0]='CHO'
        self.all_results['mol_class'][(self.all_results['C']>0) & (self.all_results['N']>0.5)]='CHON'
        #self.all_results['mol_class'][(self.all_results['C']>0) & (self.all_results['S']>0.5)]='CHOS'
        #self.all_results['mol_class'][(self.all_results['C']>0) & (self.all_results['P']>0.5)]='CHOP'
        #self.all_results['mol_class'][(self.all_results['C']>0) & (self.all_results['Na']>0.5)]='CHONa'
        #self.all_results['mol_class'][(self.all_results['C']>0) & (self.all_results['S']>0.5) & (self.all_results['N']>0.5)]='CHONS'
        #self.all_results['mol_class'][(self.all_results['C']>0) & (self.all_results['P']>0.5) & (self.all_results['N']>0.5)]='CHONP'
        #self.all_results['mol_class'][(self.all_results['C']>0) & (self.all_results['Na']>0.5) & (self.all_results['N']>0.5)]='CHONNa'
        #self.all_results['mol_class'][(self.all_results['C']>0) & (self.all_results['P']>0.5) & (self.all_results['Na']>0.5) & (self.all_results['N']>0.5)]='CHONPNa'


        results= self.all_results[self.all_results['Is Isotopologue']==0]
        results['O/C']=results['O']/results['C']
        results['H/C']=results['H']/results['C']
        results['N/C']=results['N']/results['C']


        colors = {'CHO':'red', 'CHON':'green', 'CHOS':'blue', 'CHONS':'yellow', 'CHONP':'black', 'CHONNa':'cyan','CHONPNa':'pink','CHONa':'aqua','CHOP':'gray'}

        grouped=results.groupby('mol_class')
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1,figsize=(8,12))
        ax1 = sns.violinplot(x="Time", y="O/C", hue='mol_class', data=results, ax=ax1)
        ax1.set(xlabel='Time (min)')
        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        ax2 = sns.scatterplot(x='m/z',y='m/z Error (ppm)',hue='mol_class',data=results,ax=ax2)
        ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        #for key, group in grouped:
        #    group.plot(ax=ax,kind='scatter',x='m/z',y='m/z Error (ppm)',color=colors[key],label=key)

        for key, group in grouped:
            group.plot(ax=ax3,kind='scatter',x='m/z',y='Resolving Power',color=colors[key],label=key)
        ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        assign_summary=[]
        for time in self.all_results['Time'].unique():
            current={}
            current['Time']=time
            for mol_class in sorted(self.all_results['mol_class'].unique()):
                current[mol_class]=len(self.all_results[(self.all_results['mol_class']==mol_class) & (self.all_results['Time']==time)])
            assign_summary.append(current)
        df=pd.DataFrame(assign_summary)

        ax4 = df.plot.bar(x='Time',y=df.columns[1:],stacked=True,ylabel='Peaks',ax=ax4)
        ax4.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        fig.savefig('%sall_results.pdf' %(self.savedir), bbox_to_inches='tight')

        plt.show()
        

    def determine_unique_features(self,blankfile,blnk_thresh):
        self._unique_run = True
        #Create a list of all unique features and describe their intensity. 
        print('total # results: %s' %len(self.all_results))
        #define a list of unique features (time, formula) with 'areas' determined for each sample. There may be a slight bug that causes the unique list to grow...
        uniquelist=[]
        for time in self.all_results.Time.unique():
            current=self.all_results[self.all_results.Time==time]
            current=current.sort_values(by=['m/z Error (ppm)'],ascending=True)
            currentunique=current.drop_duplicates(subset=['Molecular Formula'])
            currentunique=currentunique[currentunique['C']>1]
            currentunique=currentunique.set_index(['Molecular Formula'],drop=False)
            for file in self.all_results['file'].unique():
                current_file=current[current['file']==file].drop_duplicates(subset=['Molecular Formula'])
                current_file=current_file.rename(columns={'Peak Height':file})
                current_file=current_file.set_index(['Molecular Formula'],drop=False)
                currentunique=currentunique.join(current_file[file])
            uniquelist.append(currentunique)

        self.unique_results=pd.concat(uniquelist,ignore_index=True)
        self.unique_results['N/C']=self.unique_results['N']/self.unique_results['C']
        self.unique_results['blank']=self.unique_results[blankfile]/self.unique_results['Peak Height']
        self.unique_results=self.unique_results[self.unique_results['blank']<blnk_thresh]

        print('# unique results: %s' %len(self.unique_results))

        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1,figsize = (8,12))
        sns.violinplot(x="Time", y="O/C", hue='mol_class', data=self.unique_results, ax=ax1)
        sns.violinplot(x="Time", y="N/C", hue='mol_class', data=self.unique_results, ax=ax2)
        ax1.set(xlabel='Time (min)')
        ax2.set(xlabel='Time (min)')
        sns.scatterplot(x='m/z',y='m/z Error (ppm)',hue='N',data=self.unique_results,ax=ax3,s=10)
        sns.scatterplot(x='m/z',y='Resolving Power',hue='mol_class',data=self.unique_results,ax=ax4)
        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax4.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        fig.savefig(self.savedir + 'unique_results.pdf', bbox_to_inches='tight')
        plt.tight_layout()
        plt.show()


    def run_cluster_analysis(self, sort_order):

        #Generate clusters of ms features across depth.
        from sklearn.cluster import AgglomerativeClustering
        import seaborn as sns
        from matplotlib.colors import LogNorm, Normalize

        if self._unique_run is False:
            
            self.assess_all_results()
            self.determine_unique_features()


        only_samples = {k: self.complete_results[k] for k in self.complete_results.keys() if   ('stdmix' not in k) and ('blank' not in k)}
        self.complete_results_df=pd.concat(only_samples.values())

        #clustermethod='average'
        clustermethod='ward'

        abundances=self.unique_results[self.complete_results_df['file'].unique()].fillna(0)

        df=abundances.mean(axis=1)
        df_std=abundances.std(axis=1)

        p_list=[]
        for ind in abundances.index:
            p=max(abs(abundances.loc[ind]-df[ind]))/df_std[ind]
            if len(abundances.loc[ind][abundances.loc[ind]>1])<2:
                p=0
            p_list.append(p)

        #abundances['p']=p_list
        #abundances=abundances[abundances['p']>0.1]
        #abundances=abundances.sort_values(by='p',ascending=False)

        self.unique_results['p']=p_list

        #plt.hist(p_list,bins=100,range=[min(p_list),max(p_list)])

        results_clustered=self.unique_results[self.unique_results['p']>0.1]
        norm_abundances=results_clustered[self.complete_results_df['file'].unique()].fillna(0)
        norm_abundances=norm_abundances.div(norm_abundances.max(axis=1),axis=0)


        cluster = AgglomerativeClustering(n_clusters=3,affinity='euclidean',linkage=clustermethod)
        cluster.fit_predict(norm_abundances)

        results_clustered['cluster']=cluster.labels_

        #results_clustered.fillna(0).to_csv(file_location+'clustered_results.csv')

        clusterplot=norm_abundances
        #clusterplot=clusterplot.drop(['p'],axis=1)
        clusterplot=clusterplot.transpose()
        clusterplot['sort_order']=sort_order
        clusterplot=clusterplot.sort_values(by='sort_order')
        clusterplot=clusterplot.drop(['sort_order'],axis=1)

       # fig, ax = plt.subplots()

        p1 = sns.clustermap(clusterplot,row_cluster=False,cmap='mako',method=clustermethod)

        fig = plt.gcf()
        fig.savefig(self.savedir + 'clustering_results.pdf')

        plt.show()




