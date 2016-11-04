#!/usr/bin/env python

import numpy as np
from Fluence import UHE_fluence
from Inference import UHE_emcee
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


#parms = [norm, spectral_index, E_max, f_He, f_N, f_Si, f_Fe, source_index, z1, z2, z3, beta]
#parms = [1.33e-62,-2.2,20.0,1.11e-04,1.10e-04,1.08e-04,9.92e-02,3.0,1.87,2.66,3.01, 100.]
initial_parms = [1.20557037e-62,  -2.24416348e+00,   2.00390226e+01,
                 1.05732808e-04,   7.60038773e-06,   1.06653878e-04,
                 4.13123141e-01,   8.06863919e-01,   1.69326177e+00,
                 1.69326213e+00,   1.69326225e+00,   3.74237024e+01]

UHE_emcee(fcalc, # an instance of the UHE_fluence calculator class
          data=['Auger'], # experimental constraints
          initial_parm_vals = initial_parms,
          ndim = 12,
          nwalkers = 120, # ndim * 10
          niterations = 100000, # number of steps in mcmc sampler
          interval = 100,
          out_tag = out_tag,
          out_dir = out_dir
          )

