#!/usr/bin/env python

import numpy as np
from Fluence import UHE_fluence
import Auger_Data
from scipy.optimize import minimize
from Likelihood import UHE_likelihood
import os

crane_dir = '/home/romerowo/CRANE'
out_tag = 'run01'
out_dir = '/halo_nobackup/eva/romerowo/crane_inference_outputs/20161103'
inputs = eval(open("%s/inputs.txt"%crane_dir).read()) 

print 'Loading Lookup Table'
f = np.load(inputs['Lookup_Table'])
input_Z = f['input_Z']
input_A = f['input_A']
input_log10_E = f['input_log10_E']
input_log10_z = f['input_log10_z']
observed_log10_E = f['observed_log10_E']
observed_A = f['observed_A']
uhe_LUT = f['LUT']
fcalc = UHE_fluence(input_Z, input_A, input_log10_E, input_log10_z, observed_A, observed_log10_E, uhe_LUT)

def sample_and_minimize(N=10, lp_cut = -10000.):
  count = 0
  parms = []
  ll_vals = []
  logposterior = UHE_likelihood(fcalc, ['Auger'])
  parms_array = []
  while count < N :
    norm = 1.e-62 # this will get re-normalized before minimization
    spectral_index = np.random.uniform(-10.,0.)
    E_max = np.random.uniform(17., 25.)
    f_vals = np.random.uniform(0.,1., 5) 
    f_p, f_He, f_N, f_Si, f_Fe = f_vals / np.sum(f_vals)
    source_index = np.random.uniform(-25., 25.)
    z1 = np.random.uniform(1.e-4, 10.)
    z2 = np.random.uniform(z1, 10.)
    z3 = np.random.uniform(z2, 10.)
    beta = np.exp(np.random.uniform(-100., 100.))
    parms = [norm, spectral_index, E_max, f_He, f_N, f_Si, f_Fe, source_index, z1, z2, z3, beta]

    #print parms

    fluences = fcalc.fluence_model(*parms)
    nuclear_fluence = np.sum(fluences[1:], axis=0)
    N_model = Auger_Data.Counts_Model(nuclear_fluence[fcalc.observed_log10_E>18.8999])

    N_data = Auger_Data.Counts[Auger_Data.log10_Energy_low_edges > 18.89999]

    if( N_model[0] < 1.e-100 ):
        continue

    #print N_data[0]/N_model[0], np.sum(N_data)/np.sum(N_model)

    parms[0] *= N_data[0]/N_model[0]
    lp_init = logposterior(parms)

    #if( not np.isfinite(lp_init)):
    #    continue
    if( lp_init < lp_cut):
        continue
    
    parm_str = ''
    for k in range(0,len(parms)):
        parm_str += '%+1.2e\t'%parms[k]

    print 'lp_init %1.2f'%lp_init, parm_str

    # Check
    fluences = fcalc.fluence_model(*parms)
    nuclear_fluence = np.sum(fluences[1:], axis=0)
    N_model = Auger_Data.Counts_Model(nuclear_fluence[fcalc.observed_log10_E>18.8999])
    N_data = Auger_Data.Counts[Auger_Data.log10_Energy_low_edges > 18.89999]

    #print N_data[0]/N_model[0], np.sum(N_data)/np.sum(N_model)
    #print 'good'
    #print ''


    def minus_lnprob(_parms):
        LL = logposterior(_parms)
        #'''
        #print_str = '%1.2e '%LL
        #for k in range(0,len(_parms)):
        #    print_str += '%+1.1e '%_parms[k]
        # print the values occasionally
        #if(np.random.randint(0,100)==0): print print_str
        #'''
        return -1.*LL
    #print 'lp %1.2f '% lp_init

    res = minimize(minus_lnprob, parms, method='nelder-mead', options={'xtol': 1e-3, 'disp': True})
    parms = res.x
    lp_mini = logposterior(parms)

    parm_str = ''
    for k in range(0,len(parms)):
        parm_str += '%+1.2e\t'%parms[k]

    print 'lp_mini %1.2f'%lp_mini, parm_str
    print ''
    if( not res.success and lp_mini<-30.):
        print 'rejected'
        continue
    parms_array.append(parms)
    ll_vals.append(lp_mini)
    count += 1

    # save results
    if(count>0 and count%10 == 0):
        print 'Saving', count
        np.savez('/halo_nobackup/eva/romerowo/crane_inference_outputs/20161104/min_parms.npz',  
                 ll_array = ll_vals, 
                 parms_array = np.array(parms_array))


  np.savez('/halo_nobackup/eva/romerowo/crane_inference_outputs/20161104/min_parms.npz',  
           ll_array = ll_vals, 
           parms_array = np.array(parms_array))

sample_and_minimize(N=600, lp_cut = -50.)
