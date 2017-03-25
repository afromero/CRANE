import numpy as np
#import emcee
#import uhe_data
#from scipy.stats import poisson
#from scipy.special import factorial

import Auger_Data
from Xmax import Xmax_modeler

class UHE_likelihood(object):
    # initialize by passing the UHE_fluence calculator
    # the data input is a list that specifies which data to use ['Auger', 'Auger_Xmax', 'IceCube', 'Auger_Neutrino', etc...]
    def __init__(self, UHE_fluence, data, 
        log10_E_sys = 0.0569,  # log10(1+0.14) corresponding to 14% systematic reported by Auger 
        norm_lower = 0., norm_upper = np.inf,
        spectral_index_lower = -np.inf, spectral_index_upper = 0.,
        E_max_lower = 10., E_max_upper = 25.,
        f_p_lower   = 0.,  f_p_upper   = 1.,
        f_He_lower  = 0.,  f_He_upper  = 1.,
        f_N_lower   = 0.,  f_N_upper   = 1.,
        f_Si_lower  = 0.,  f_Si_upper  = 1.,
        f_Fe_lower  = 0.,  f_Fe_upper  = 1.,
        source_index_lower = -np.inf, source_index_upper = np.inf,
        z1_lower = 0.,     z1_upper   = 10.,
        z2_lower = 0.,     z2_upper   = 10.,
        z3_lower = 0.,     z3_upper   = 10.,
	log10_E_shift_lower = -np.inf, log10_E_shift_upper=+np.inf,
        uX_lower = 0,      uX_upper   = 1.):
        #self.tmax = 10.0*np.pi
        #self.constant = np.log(1.0/(self.tmax*self.tmax))
        # can define bounds here, can also pass data to initialize the function?
        print 'Initializing UHE_likelihood'
        self.UHE_fluence  = UHE_fluence
        self.data         = data
        self.Xmax_log10_E_lower_edges = np.array([18.9, 19.0, 19.1, 19.2, 19.3, 19.4, 19.5])
        self.Xmax_modeler = Xmax_modeler(self.Xmax_log10_E_lower_edges)
        
        self.log10_E_sys          = log10_E_sys

        self.norm_lower           = norm_lower
        self.norm_upper           = norm_upper
        self.spectral_index_lower = spectral_index_lower
        self.spectral_index_upper = spectral_index_upper
        self.E_max_lower          = E_max_lower
        self.E_max_upper          = E_max_upper
        self.f_p_lower            = f_p_lower
        self.f_p_upper            = f_p_upper
        self.f_He_lower           = f_He_lower
        self.f_He_upper           = f_He_upper
        self.f_N_lower            = f_N_lower
        self.f_N_upper            = f_N_upper
        self.f_Si_lower           = f_Si_lower
        self.f_Si_upper           = f_Si_upper
        self.f_Fe_lower           = f_Fe_lower
        self.f_Fe_upper           = f_Fe_upper
        self.source_index_lower   = source_index_lower
        self.source_index_upper   = source_index_upper
        self.z1_lower             = z1_lower
        self.z1_upper             = z1_upper
        self.z2_lower             = z2_lower
        self.z2_upper             = z2_upper
        self.z3_lower             = z3_lower
        self.z3_upper             = z3_upper
        self.log10_E_shift_lower  = log10_E_shift_lower
	self.log10_E_shift_upper  = log10_E_shift_upper
        self.uX_lower             = uX_lower
        self.uX_upper             = uX_upper

    # Set bounds and prior distributions on the parameters
    def logprior(self, _parms):
        # read the fluence model parameter array.
        # norm, spectral_index, E_max, f_He, f_N, f_Si, f_Fe, source_index, z1, z2, z3, beta = _parms (20-Mar-2017, removing beta Fermi function parameter in favor of an exponential)
        # if( (not np.isfinite(beta)) or (beta<=0.) or (beta>3.e100)):
        #     return -np.inf
        norm, spectral_index, E_max, f_He, f_N, f_Si, f_Fe, source_index, z1, z2, z3, log10_E_shift, uX = _parms
 
       ####################################
       # some pretty basic priors for now #
       ####################################
    
        # the norm range
        pass_cond = norm>self.norm_lower
        pass_cond = np.logical_and(norm < self.norm_upper, pass_cond)
        #print ''
        #print 'pass_cond norm', pass_cond
        
        # The source spectral index range
        pass_cond = np.logical_and(spectral_index > self.spectral_index_lower, pass_cond)
        pass_cond = np.logical_and(spectral_index < self.spectral_index_upper, pass_cond)
        #print 'pass_cond spec_ind', pass_cond

        # The E_max range
        pass_cond = np.logical_and(E_max > self.E_max_lower, pass_cond)
        pass_cond = np.logical_and(E_max < self.E_max_upper, pass_cond)
        #print 'pass_cond E_max', pass_cond

        # The source compositon fractions must be positive and their sum is smaller than unity.
        f_p = 1. - f_He - f_N - f_Si - f_Fe
        pass_cond = np.logical_and(f_p  >= self.f_p_lower, pass_cond)
        pass_cond = np.logical_and(f_p  <= self.f_p_upper, pass_cond)
        #print 'pass_cond f_p', pass_cond

        pass_cond = np.logical_and(f_He >= self.f_He_lower, pass_cond)
        pass_cond = np.logical_and(f_He <= self.f_He_upper, pass_cond)
        #print 'pass_cond f_He', pass_cond

        pass_cond = np.logical_and(f_N >= self.f_N_lower, pass_cond)
        pass_cond = np.logical_and(f_N <= self.f_N_upper, pass_cond)
        #print 'pass_cond f_N', pass_cond

        pass_cond = np.logical_and(f_Si >= self.f_Si_lower, pass_cond)
        pass_cond = np.logical_and(f_Si <= self.f_Si_upper, pass_cond)
        #print 'pass_cond f_Si', pass_cond

        pass_cond = np.logical_and(f_Fe >= self.f_Fe_lower, pass_cond)
        pass_cond = np.logical_and(f_Fe <= self.f_Fe_upper, pass_cond)
        #print 'pass_cond f_Fe', pass_cond

        # the sum of nuclear fractions heavier than proton cannot exceed 1. 
        pass_cond = np.logical_and(f_He + f_N + f_Si + f_Fe <= 1., pass_cond)
        #print 'pass_cond f_Sum', pass_cond

        # The source distribution parameters
        pass_cond = np.logical_and(source_index > self.source_index_lower, pass_cond)
        pass_cond = np.logical_and(source_index < self.source_index_upper, pass_cond)
        #print 'pass_cond source_index', pass_cond

        pass_cond = np.logical_and(z1 > self.z1_lower, pass_cond)
        pass_cond = np.logical_and(z1 < self.z1_upper, pass_cond)
        #print 'pass_cond z1', pass_cond

        pass_cond = np.logical_and(z2 > self.z2_lower, pass_cond)
        pass_cond = np.logical_and(z2 < self.z2_upper, pass_cond)
        #print 'pass_cond z2', pass_cond

        pass_cond = np.logical_and(z3 > self.z3_lower, pass_cond)
        pass_cond = np.logical_and(z3 < self.z3_upper, pass_cond)
        #print 'pass_cond z3', pass_cond

        # The source distribution inflection points must increase monotonically and be within the bounds of the Yields calculation.
        pass_cond = np.logical_and(z2>z1, pass_cond)
        pass_cond = np.logical_and(z3>z2, pass_cond)
        #print 'pass_cond monotonic', pass_cond
      
        pass_cond = np.logical_and(log10_E_shift > self.log10_E_shift_lower, pass_cond)
        pass_cond = np.logical_and(log10_E_shift < self.log10_E_shift_upper, pass_cond)

        pass_cond = np.logical_and(uX >= self.uX_lower, pass_cond)
        pass_cond = np.logical_and(uX <= self.uX_upper, pass_cond)

        # At some point we may want to include tighter bounds on metallicity or other restriction

        # If all conidtions pass, the priors have probability equalt to unity.
        # If a condition fails, the probability is zero.
        if(pass_cond==True):
            return -0.5*log10_E_shift**2/self.log10_E_sys**2
        return -np.inf
        

    def loglhood(self,_parms):
        fluence_parms = _parms[0:11]
	log10_E_shift = _parms[11]
        # Calculate fluence curve for these parameters
        fluences = self.UHE_fluence.fluence_model(*fluence_parms)
        nuc_fluence = np.sum(fluences[1:], axis=0)
        
        LL = 0. # initialize likelihood
        #print fluences[0,:]
        if 'Auger' in self.data:
            # model Auger counts
            N_model         = Auger_Data.Counts_Model(nuc_fluence[self.UHE_fluence.observed_log10_E>18.8999-0.2], log10_E_shift)
            # use energy of 18.69999 so as to have 2 bins of margin (needed for energy systematics)
            N_model = N_model[2:]
            N_data          = Auger_Data.Counts[Auger_Data.log10_Energy_low_edges > 18.8999]
            N_data_log_fac  = Auger_Data.log_factorial[Auger_Data.log10_Energy_low_edges > 18.8999]
	    
            #print 'N_model, N_data lengths = ', len(N_model), len(N_data)
            log_poisson = np.zeros(len(N_model))
            log_poisson[N_model>0.] = N_data[N_model>0.]*np.log(N_model[N_model>0.])  - N_model[N_model>0.] - N_data_log_fac[N_model>0.]

            # handle the case where zero counts are expected and zero are obtained.
            # This has probability = 1.
            log_poisson[np.logical_and(N_model==0, N_data==0)] = 0.

            # handle the case where zero counts are expected and more than zero counts occured
            # such a model should definitely be rejected (with the caveat that energy bin spill-over could do this)
            log_poisson[np.logical_and(N_model==0,N_data!=0)] = -np.infty
          
            LL += np.sum(log_poisson)

        if 'Auger_Xmax' in self.data:
            fluence_energy_bounds = np.logical_and(self.UHE_fluence.observed_log10_E > np.min(self.Xmax_log10_E_lower_edges) - 1.e-3, self.UHE_fluence.observed_log10_E < np.max(self.Xmax_log10_E_lower_edges) + 1.e-3) 
            data_energy_bounds    = np.logical_and(Auger_Data.X_max_log10_E_low_edge > np.min(self.Xmax_log10_E_lower_edges) - 1.e-3, Auger_Data.X_max_log10_E_low_edge< np.max(self.Xmax_log10_E_lower_edges) + 1.e-3) 
            f_A_array = fluences[1:, fluence_energy_bounds]/nuc_fluence[fluence_energy_bounds]
            f_A_array = np.nan_to_num(f_A_array) # if nuc_fluence is zero then f_A_array is zero
            
            #model_Mean, model_RMS = self.Xmax_modeler.getMeanRMS(f_A_array)  # best not to precalculate Xmax_modeler if we want to include energy systematics

            # we want to get the Xmax values at the true energy but compare them to the observed energy (by the detector)
            # observed_energy = real_energy + log10_E_shift, so the real_energy = observed_energy - log10_E_shift

            xmm = Xmax_modeler(self.Xmax_log10_E_lower_edges - log10_E_shift) 
            model_Mean, model_RMS = xmm.getMeanRMS(f_A_array) 

            #model_Mean = np.mean(model_Mean, axis=0)
            #model_RMS  = np.mean(model_RMS,  axis=0)
            uX = _parms[12]
            model_Mean  = uX*model_Mean[0] + (1.-uX)*model_Mean[1]
            model_RMS   = uX*model_RMS[0]  + (1.-uX)*model_RMS[1]
            

            LL += np.sum( -0.5*((model_Mean - Auger_Data.X_max_Mean[data_energy_bounds])/Auger_Data.X_max_Mean_err[data_energy_bounds] )**2 )
            LL += np.sum( -0.5*((model_RMS  - Auger_Data.X_max_RMS[data_energy_bounds])/Auger_Data.X_max_RMS_err[data_energy_bounds] )**2 )

        return LL
    def __call__(self, _parms):
        lprior = self.logprior(_parms)
        # don't waste time calculating likelihoods for parameters that are not allowed.
        if( np.isinf(lprior)):
            return lprior

        # return the prior plus the calculated log likelihood
        if(np.isnan(lprior)):
            print 'lprior is nan'
        if(np.isnan(self.loglhood(_parms))):
            print 'loglhood is nan'
        return  lprior + self.loglhood(_parms)


