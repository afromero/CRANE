import numpy as np
#import emcee
#import uhe_data
#from scipy.stats import poisson
#from scipy.special import factorial

import Auger_Data

class UHE_likelihood(object):
    # initialize by passing the UHE_fluence calculator
    # the data input is a list that specifies which data to use ['Auger', 'Auger_Xmax', 'IceCube', 'Auger_Neutrino', etc...]
    def __init__(self, UHE_fluence, data):
        #self.tmax = 10.0*np.pi
        #self.constant = np.log(1.0/(self.tmax*self.tmax))
        # can define bounds here, can also pass data to initialize the function?
        print 'Initializing UHE_likelihood'
        self.UHE_fluence = UHE_fluence
        self.data = data
        
    # Set bounds and prior distributions on the parameters
    def logprior(self, _parms):
        # read the fluence model parameter array.
        norm, spectral_index, E_max, f_He, f_N, f_Si, f_Fe, source_index, z1, z2, z3, beta = _parms
        if( (not np.isfinite(beta)) or (beta<=0.) or (beta>3.e100)):
            return -np.inf
 
       ####################################
       # some pretty basic priors for now #
       ####################################
    
        # the norm is positive
        pass_cond = norm>0.
        
        # The source spectral index is negative
        pass_cond = np.logical_and(spectral_index < 0., pass_cond)

        # The maximum energy of acceleration is within the bounds provided by the Yields calcluation
        pass_cond = np.logical_and(E_max < 25., pass_cond)

        # The source compositon fractions must be positive and their sum is smaller than unity.
        pass_cond = np.logical_and(f_He >= 0., pass_cond)
        pass_cond = np.logical_and(f_N  >= 0., pass_cond)
        pass_cond = np.logical_and(f_Si >= 0., pass_cond)
        pass_cond = np.logical_and(f_Fe >= 0., pass_cond)
        pass_cond = np.logical_and(f_He + f_N + f_Si + f_Fe <= 1., pass_cond)

        # The source distribution inflection points must increase monotonically and be within the bounds of the Yields calculation.
        pass_cond = np.logical_and(z1>1.e-4, pass_cond)
        pass_cond = np.logical_and(z2>z1, pass_cond)
        pass_cond = np.logical_and(z3>z2, pass_cond)
        pass_cond = np.logical_and(z3<=10., pass_cond)
      
        # At some point we may want to include tighter bounds on metallicity or other restriction

        # If all conidtions pass, the priors have probability equalt to unity.
        # If a condition fails, the probability is zero.
        if(pass_cond==True):
            return 0.
        return -np.inf
        

    def loglhood(self,_parms):
        # Calculate fluence curve for these parameters
        fluences = self.UHE_fluence.fluence_model(*_parms)
        
        LL = 0. # initialize likelihood
        #print fluences[0,:]
        if 'Auger' in self.data:
            # model Auger counts
            nuc_fluence = np.sum(fluences[1:], axis=0)
            N_model         = Auger_Data.Counts_Model(nuc_fluence[self.UHE_fluence.observed_log10_E>18.8999])
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
        return LL
        
    def __call__(self, _parms):
        lprior = self.logprior(_parms)
        # don't waste time calculating likelihoods for parameters that are not allowed.
        if( np.isinf(lprior)):
            return lprior

        # return the prior plus the calculated log likelihood
        return  lprior + self.loglhood(_parms)


