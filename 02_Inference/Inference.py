from Likelihood import UHE_likelihood
import numpy as np
import emcee
import pickle
import sys

def dither_parms(parms):
    # Dither the initial positions
    parms       = parms*(1. + 1.e-8*np.random.randn(len(parms)))
    parms[3:7] += np.random.uniform(0.,1.e-8,4)
    return parms

def init_pos_from_npz_file(init_npz_filename, pos):
    f = np.load(init_npz_filename)
    #print f.keys()
    ll_array    = f['ll_array']
    parms_array = f['parms_array']
    f.close()
    for k in range(0,len(pos)):
        parms = parms_array[np.random.randint(0,len(parms_array))]
        #print ''
        #print parms
        pos[k] = dither_parms(parms)
        #print pos[k]
    #print pos.shape
    return pos
    

def UHE_emcee(  UHE_fluence, # an instance of the UHE_fluence calculator class
                data=['Auger', 'Auger_Xmax', 'ANITA', 'IceCube', 'ARA', 'EVA'], # experimental constraints
                initial_parm_vals = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., 0.],
                init_npz_filename = None, # if a file name is given, then it will read it.
                ndim = 12,
                nwalkers = 120, # ndim * 10
                niterations = 1000, # number of steps in mcmc sampler
                interval = 10, # the interval over which we save the data.
                out_tag = '',
                out_dir = './',
                log10_E_sys = 0.0569,  # log10(1+0.14) corresponding to 14% systematic reported by Auger 
                norm_lower = 0., norm_upper = np.inf,
                spectral_index_lower = -np.inf, spectral_index_upper = +np.inf,
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
    
    # initialize the likelihood function
    logposterior = UHE_likelihood(UHE_fluence, data,
                              norm_lower = norm_lower, 
                              norm_upper = norm_upper,
                              spectral_index_lower = spectral_index_lower, 
                              spectral_index_upper = spectral_index_upper,
                              E_max_lower = E_max_lower, 
                              E_max_upper = E_max_upper,
                              f_p_lower   = f_p_lower ,  
                              f_p_upper   = f_p_upper,
                              f_He_lower  = f_He_lower,  
                              f_He_upper  = f_He_upper,
                              f_N_lower   = f_N_lower,  
                              f_N_upper   = f_N_upper,
                              f_Si_lower  = f_Si_lower,  
                              f_Si_upper  = f_Si_upper,
                              f_Fe_lower  = f_Fe_lower,  
                              f_Fe_upper  = f_Fe_upper,
                              source_index_lower = source_index_lower, 
                              source_index_upper = source_index_upper,
                              z1_lower = z1_lower,     
                              z1_upper = z1_upper,
                              z2_lower = z2_lower,     
                              z2_upper = z2_upper,
                              z3_lower = z3_lower,     
                              z3_upper = z3_upper,
                              log10_E_sys = log10_E_sys,  
                              log10_E_shift_lower = log10_E_shift_lower,
                              log10_E_shift_upper = log10_E_shift_upper,
                              uX_lower = uX_lower,     
                              uX_upper = uX_upper)

    print 'Initial Likelihood', logposterior(initial_parm_vals)

    # Set up the emcee Ensemble sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, logposterior, threads=1)
    
    # Dither the initial positions
    pos = np.array( [initial_parm_vals*(1. + 1.e-2*np.random.randn(ndim)) for i in range(nwalkers)] )
    # Add a little bit of composition
    for k in range(0,len(pos)):
        pos[k][3:7] += np.random.uniform(0.,0.01,4)

    if( init_npz_filename ):
        #print 'init_npz_filename not None!'
        pos = init_pos_from_npz_file(init_npz_filename, pos)
        #exit()

    #for k in range(0,len(pos)):
    #    np.set_printoptions(precision=3)
    #    print logposterior(pos[k]), pos[k] 

    #exit()

    sys.stdout.flush()
    # Set up an interval of iterations over which to output the chains.
    N_loops = niterations / interval
    print interval, N_loops, N_loops*interval, niterations
    for N in range(0,N_loops):
        print 'N=', N
        sys.stdout.flush()
        # SAVING THE 
        #pickle.dump(sampler, open("uhe_emcee_sampler.p", "wb"))
        pos, prob, state = sampler.run_mcmc(pos, interval)
        out_chain  = sampler.chain[:,N*interval:(N+1)*interval,:]
        out_lnprob = sampler.lnprobability[:,N*interval:(N+1)*interval]

        #pickle.dump(sampler.chain, open("uhe_emcee_chain_%d.p"%N, "wb"))
        #pickle.dump(sampler.lnprobability, open("uhe_emcee_lnprobability_%d.p"%N, "wb"))
        print '\tsampler.chain.shape', sampler.chain.shape, sampler.lnprobability.shape

        print '\tout_chain.shape', out_chain.shape, out_lnprob.shape
        pickle.dump(out_chain, open("%s/uhe_emcee_chain_%s_%d.p"%(out_dir,out_tag,N), "wb"))
        pickle.dump(out_lnprob, open("%s/uhe_emcee_lnprob_%s_%d.p"%(out_dir,out_tag,N), "wb"))
        print '\t%1.2e'%(prob[0])
        print '\t%1.2e\t%1.1f\t%1.1f\t%1.2e\t%1.2e\t%1.2e\t%1.2e\t%1.1f\t%1.2f\t%1.2f\t%1.2f\t%1.2f'%(pos[0][0], pos[0][1], pos[0][2], pos[0][3], pos[0][4], pos[0][5], pos[0][6], pos[0][7], pos[0][8], pos[0][9], pos[0][10], pos[0][11])
        # Print out the mean acceptance fraction. In general, acceptance_fraction
        print "\tMean acceptance fraction:", np.mean(sampler.acceptance_fraction)
        # Estimate the integrated autocorrelation time for the time series in each
        # parameter.
        print "\tAutocorrelation time:", sampler.get_autocorr_time()
        sys.stdout.flush()

