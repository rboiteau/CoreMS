
    
if __name__ == '__main__':
    startdt = datetime.now()
    start = time.time()  #for duration
   
    drive_dir = '/Volumes/Samsung_T5/ESI-ICP-MS/seawater-iodine/'
    svdir=drive_dir
    
    mzref = "/Users/christiandewey/CoreMS/db/Hawkes_neg.ref"
    
    heteroAtom = '127I'

    offset = -36 # ICPMS time (s) + offset (s) = ESI time (s)




    esifile_name = '220822_CTD27_600m2.raw'
    data_dir = drive_dir
    esifile = data_dir+esifile_name
    esiparser = rawFileReader.ImportMassSpectraThermoMSFileReader(esifile)
    icpmsfile = data_dir + 'CTD27_600m.csv'
    results_fname='230301_CWD_220822_CTD27_600m2_001.csv'

      
    '''### plot offset ICPMS data and select peak for matching 
                fig, ax = plt.subplots()
                global coords
                coords = []
                cid = fig.canvas.mpl_connect('button_press_event', mouse_event)
                ax = plotICPMS(icpmsfile,['127I', '59Co'], offset, ax)
                ax.set_xlim(0,10)
                plt.show()
            
                
                ### subset ICP data 
                trange = [coords[0][0], coords[1][0]]'''

    trange = [44,47]  



    icpsub = subset_icpdata(icp_data_file=icpmsfile, heteroAtom='127I', timerange=trange, offset = offset)
  
    
    ### interpolate ICP data to match ESI data time points
    interpolated_ICP = interpolate(esi_parser=esiparser,icpsub=icpsub,heteroAtom='127I',timerange=trange)

    
    ### run formula assignment 
    setAssingmentParams()
    print(os.getcwd())
    assignments = assign_formula(esiparser,trange,mzref,cal_ppm_threshold=(-10,10),charge=1)
    assignments.to_csv(svdir+results_fname)
    all_results = pd.read_csv(svdir+results_fname)


    ### get EICS
    EICs, avMS = get_eics(esi_parser=esiparser,assignments = assignments, timerange=trange)



    with open(svdir + 'eics_'+results_fname + '.p', 'wb') as fp:
                                pickle.dump(EICs, fp, protocol=pickle.HIGHEST_PROTOCOL)
            
    with open(svdir + 'eics_'+results_fname + '.p', 'rb') as fp:
        EICs = pickle.load(fp)
    ### run correlation and filter assignment results
    corr_assignments, corr_assignments_with_hetero = correlate(icp_interp=interpolated_ICP,EICdic=EICs,heteroAtom='127I',assignments=assignments,timerange=trange,threshold = 0.2) 

    corr_assignments.to_csv(svdir+'corr_'+results_fname)
    
    for mz in corr_assignments['m/z']:
        fig, ax = plt.subplots()
        eic_df=pd.DataFrame({'EIC':EICs[mz].eic,'Time':EICs[mz].time})
        eic_sub=eic_df[eic_df['Time'].between(trange[0],trange[1])]
        #def plot_EIC_ICPMS(eic,mz,icp_data,element, trange, ax = None):
        plot_EIC_ICPMS(eic = eic_sub,mz = mz,icp_data = interpolated_ICP, element='127I', trange=trange, ax=ax)
        mf = corr_assignments[corr_assignments['m/z'] == mz]['Molecular Formula'].iloc[0]

        if pd.isna(mf):
            mf = 'Unassigned'
            mz_err = ''
        else:
            mz_err = 'm/z Error (ppm): %.3f' %(corr_assignments[corr_assignments['m/z'] == mz]['m/z Error (ppm)'].iloc[0])
        ax.text(0.05,0.95, mf, transform=ax.transAxes)
        ax.text(0.05,0.90, mz_err, transform=ax.transAxes)
        plt.savefig(svdir + results_fname.split('.')[0] + '_' + mf + '_' +'%.4f.pdf' %(mz))

    elapsed_time_sec = (time.time() - start) 

    if elapsed_time_sec > 3600:
        elapsed_time = elapsed_time_sec/3600
        unit = 'hr'
    else:
        elapsed_time = elapsed_time_sec / 60
        unit = 'min'

    enddt = datetime.now()
    startdt_str = startdt.strftime("%d-%b-%Y %H:%M:%S")
    enddt_str = enddt.strftime("%d-%b-%Y %H:%M:%S")

    print('\nAssignment took %.2f %s to complete' %(elapsed_time, unit))
    print('Started ' + startdt_str )
    print('Finished ' + enddt_str )














