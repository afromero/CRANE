# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 16:33:17 2016
@author: romerowo
"""

from pylab import *

# Constants
speed_of_light = 2.9979e10  # cm /sec
#H0 = 1./12.e3  # Myr^-1
H0 = 1./4.56e17  # s^-1
Omega_M = 0.3
Omega_L = 0.7
erg2eV = 6.242e+11 # is this right?
P0 = 4.5e41 # erg/ yr
Mpc2cm  = 3.085e+24

class UHE_fluence:
    def __init__(self, input_Z, input_A, input_log10_E, input_log10_z, observed_A, observed_log10_E, uhe_LUT):
        print 'Initializing UHE_fluence'
        # initiate the look-up table arrays in the class
        self.input_Z = input_Z
        self.input_A = input_A
        self.input_log10_E = input_log10_E
        self.input_log10_z = input_log10_z
        self.observed_A = observed_A
        self.observed_log10_E = observed_log10_E
        self.uhe_LUT = uhe_LUT

        # define some auxiliary arrays to avoid calculating them every iteration
        self.z = np.array(10**input_log10_z)
        self.diff_array = self.differential_array()


    def differential_array(self): 
        # this calculates the differentials for the model array integration
        # only needs to be computed once since it does not change
        # the units of the differential array are in cm * eV 
        dA = np.ones(len(self.input_Z))
        dE = np.log(10.) * 10**self.input_log10_E * diff(self.input_log10_E)[0]
        dz = np.log(10.) * self.z * diff(self.input_log10_z)[0]
        cosmo = speed_of_light / H0 / (1.+ self.z) / np.sqrt(Omega_M*(1+self.z)**3 + Omega_L)
        diff_array = np.einsum('i,j,k->ijk', dA, dE, cosmo*dz)
        return diff_array

    def source_energy_spectrum(self,spectral_index=-2, log10_E_ref=18.):
        # this function returns the differential source energy spectrum for the integrand.
        src_flux = (10.**(self.input_log10_E - log10_E_ref) )**spectral_index
        
        # not including an E_min cut.
        # src_flux[self.input_log10_E < log10_E_min] *= 0. 

        # There are some concerns here. 
        # 1. The normalization of the source flux will depend on the energy cutoff, which also depends on the charge.
        # 2. This affects the normalization of the flux, which does have a physical meaning. 
        # 3. Some more thought is needed to figure out the best way to normalize. 
        #    It may be best done after the E_max cutoff and do it for each charge, then apply the relative abundance.
        #
        return src_flux
    
    def source_distribution(self,source_index=3., z1=1.9, z2=2.7, z3=5.):
        # initialize the array
        source_pop = np.zeros(len(self.z))

        # piecewise function with a source index power law for z between 0 and z1, 
        # flat between z1 and z2, 
        # and exponentially falling from z2 to z3
        # hard cutoff at z3
        source_pop[self.z<=z1] = (1.+self.z[self.z<=z1])**source_index
        source_pop[np.logical_and(self.z>z1,self.z<z2)] = (1.+z1)**source_index 
        source_pop[self.z>=z2] = (1.+z1)**source_index * exp((z2-self.z[self.z>=z2])/z2)      
        source_pop[self.z>=z3] = 0.  
        return source_pop
    
    def source_model_array(self,relative_abundances, source_energy_spec, source_distrib, E_max = 22.):
        # the proton fraction is whatever is left from the sum of the rest
        relative_abundances = concatenate([[1.-np.sum(relative_abundances)], relative_abundances])

        # take the outer product of the relative abundances and source energy spectrum arrays
        model_array = np.outer(relative_abundances, source_energy_spec)

        # apply charge-dependent maximum energy cutoff
        #model_array[np.outer(1./self.input_Z, 10**self.input_log10_E) > 10**E_max] *= 0.
        # arg = np.outer(1./self.input_Z, 10**(self.input_log10_E - E_max)) * beta  - beta # argument to the Fermi function (20-Mar-2017, removing it in favor of an exponential)
        arg = np.outer(1./self.input_Z, 10**(self.input_log10_E - E_max))  # argument to the exponential function (added on 20-Mar-2017)
        #model_array *= 1. / ( 1. + np.exp(arg)  ) # writing it this way produces overflow problems
        arg[arg > 709.7]  = 709.7  # cap at exponential overflow.
        arg[arg < -709.7] = -709.7 # cap at exponential underflow.
        exp_val = np.nan_to_num(np.exp(-arg))
        #print np.min(arg), np.max(arg)
        #print np.min(exp_val), np.max(exp_val)
        # model_array *= 1. / (1. + exp_val) # multiply by Fermi function (20-Mar-2017, removing it in favor of an exponential)
        model_array *= exp_val # multiply by Fermi function
        #model_array *= np.divide(exp_val , np.add(1., exp_val) ) # np.add and np.divide handles division by infinity
        
        # include the outer product with redshift evolution
        model_array = np.einsum('ij,k->ijk', model_array, source_distrib)
        
        # in the futuere we may want to include some metalicity evolution with redshift here.
        
        return model_array
    
    def fluence_calculation(self, model_array): # this returns E_{obs}*Flux
        # multiply the model_array with the differential array, 
        # then take the inner product with the yields lookup table input axes
        return np.einsum('ijk,ijklm', model_array*self.diff_array, self.uhe_LUT)*np.log(10.)
        
    def fluence_model(self,norm, 
                           spectral_index=-2., 
                           E_max=22., 
                           f_He=0., f_N=0., f_Si=0., f_Fe=0., 
                           source_index=-3., 
                           z1=1.9, z2=2.7, z3=5.,
                           # beta=2., (20-Mar-2017, removing it in favor of an exponential)
                           log10_E_ref=18.):
        e_spec   = self.source_energy_spectrum(spectral_index=spectral_index, log10_E_ref=log10_E_ref)
        src_dist = self.source_distribution(source_index=source_index, z1=z1, z2=z2, z3=z3)
        rel_ab = [ f_He, f_N, f_Si, f_Fe] # abundance of He, N, Si, Fe. The remainder is protons
        #model_array = self.source_model_array(rel_ab, e_spec, src_dist, E_max = E_max, beta = beta) # (20-Mar-2017, removing beta Fermi function parameter in favor of an exponential)
        model_array = self.source_model_array(rel_ab, e_spec, src_dist, E_max = E_max) 
        return norm*self.fluence_calculation(model_array)

