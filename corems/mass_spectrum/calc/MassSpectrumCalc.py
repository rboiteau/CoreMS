__author__ = "Yuri E. Corilo"
__date__ = "Jun 27, 2019"

from numpy import power, multiply, sqrt, multiply, array, mean
from corems.mass_spectrum.calc.NoiseCalc import NoiseThresholdCalc
from corems.mass_spectrum.calc.PeakPicking import PeakPicking

class MassSpecCalc(PeakPicking, NoiseThresholdCalc ):
    '''
    Class including numerical calculations related to mass spectrum class
    Inherited PeakPicking and NoiseThresholdCalc ensuring its methods are 
    available to the instantiated mass spectrum class object
    '''
    def percentile_assigned(self, report_error=False):
        
        assign_abun = 0
        not_assign_abun = 0
        i = 0
        j = 0
        if report_error:
            error = []
        for mspeak in self.sort_by_abundance():
            
            if mspeak.is_assigned:
                i += 1
                assign_abun += mspeak.abundance
                if report_error:
                    error.append(mspeak.best_molecular_formula_candidate.mz_error)

            else:
                j += 1
                not_assign_abun += mspeak.abundance

        total_percent = (i / (i + j)) * 100
        total_relative_abundance = (assign_abun / (not_assign_abun + assign_abun)) * 100
        if report_error:
            rms_error = sqrt(mean(array(error)**2))
            print('%i assigned peaks and %i unassigned peaks, total  = %.2f %%, relative abundance = %.2f %%, RMS error (best candidate) (ppm) = %.3f' % (i, j, total_percent, total_relative_abundance, rms_error))
            return i, j, total_percent, total_relative_abundance, rms_error

        else:
            print('%i assigned peaks and %i unassigned peaks , total  = %.2f %%, relative abundance = %.2f %%' % (i, j, total_percent, total_relative_abundance,))
            return i, j, total_percent, total_relative_abundance

    def resolving_power_calc(self, B, T):
        '''
        low pressure limits, 
        T: float 
            transient time
        B: float
            Magnetic Filed Strength (Tesla)    
        
        reference
        Marshall et al. (Mass Spectrom Rev. 1998 Jan-Feb;17(1):1-35.)
        DOI: 10.1002/(SICI)1098-2787(1998)17:1<1::AID-MAS1>3.0.CO;2-K
        
        '''
       
        self.check_mspeaks()
        return array([mspeak.resolving_power_calc(B, T) for mspeak in self.mspeaks])
        
    def _f_to_mz(self):
        ''' Ledford equation for converting frequency(Hz) to m/z, 
        
        Attributes
        ----------
        All Attributes are derivative from the MassSpecBase Class
        by calling self
        
        Returns 
        ----------
            numpy.array(float)
            m/z domain after conversion from frequency
        '''
        Aterm, Bterm, Cterm = self.Aterm, self.Bterm, self.Cterm
        # Check if the Bterm of Ledford equation scales with the ICR trap voltage or not then Bterm = Bterm*trap_voltage
        
        if Cterm == 0:
            
            if Bterm == 0:
                #uncalibrated data
                mz_domain = Aterm / self.freq_exp_profile 
                
            else:
                
                mz_domain = (Aterm / (self.freq_exp_profile)) + (Bterm / power((self.freq_exp_profile), 2))

        # @will I need you insight here, not sure what is the inverted ledford equation that Bruker refers to
        else:

            mz_domain = (Aterm / self.freq_exp_profile) + (Bterm / power(self.freq_exp_profile, 2)) + Cterm

        return mz_domain

    def _f_to_mz_bruker(self):
        ''' 
        Burker equations for converting frequency (Hz) to m/z, 
        nOmega acquisition is not yet implemented here
        
        Attributes
        ----------
        All Attributes are derivative from the MassSpecBase Class
        
        Returns 
        ----------
            numpy.array(float)
            m/z domain after conversion from frequency
        '''
        Aterm, Bterm, Cterm = self.Aterm, self.Bterm, self.Cterm
        # Check if the Bterm of Ledford equation scales with the ICR trap voltage or not then Bterm = Bterm*trap_voltage
        if Cterm == 0:
            
            if Bterm == 0:
                #uncalibrated data
                return Aterm / self.freq_exp_profile 
            
            else:
                #calc2
                return Aterm / (self.freq_exp_profile + Bterm)

        # @will I need you insight here, not sure what is the inverted ledford equation that Bruker refers to
        else:
            diff = Aterm * Aterm
            
            #this sign(diff + 4) changes on older aquistion software
            diff = diff + 4 * Cterm * (self.freq_exp_profile - Bterm)
            diff = sqrt(diff)
            diff = -Aterm+diff
            #calc3
            return (2*Cterm)/diff
            return diff/2* (self.freq_exp_profile - Bterm)

    def number_average_molecular_weight(self, profile=False):
        ''' 
        Average molecular weight calculation 
        
        Attributes
        ----------
        All Attributes are derivative from the MassSpecBase Class
        
        Returns 
        ----------
            (float)
        '''
        # mode is profile or centroid data
        if profile:
            a = multiply(self.mz_exp_profile, self.abundance_profile)
            b = self.abundance_profile
            return a.sum()/b.sum()

        else:

            return sum(self.mz_exp*self.abundance)/sum(self.abundance)
    
    def weight_average_molecular_weight(self, profile=False):
        ''' 
        Weighted Average molecular weight calculation 
        
        Attributes
        ----------
        All Attributes are derivative from the MassSpecBase Class
        
        Returns 
        ----------
            (float)
        '''
        
        # implement from MassSpectralPeaks objs

        if profile:
            a = multiply(power(self.mz_exp_profile, 2), self.abundance_profile)
            b = self.mz_exp_profile*self.abundance_profile
            return a.sum() / b.sum()

        else:
            return sum(power(self.mz_exp, 2)*self.abundance)/sum(self.mz_exp*self.abundance)
