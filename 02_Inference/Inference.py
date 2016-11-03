from Likelihood import UHE_likelihood
import numpy as np
import emcee
import pickle

def UHE_emcee(  UHE_fluence, # an instance of the UHE_fluence calculator class
                data=['Auger', 'ANITA', 'IceCube', 'ARA', 'EVA'], # experimental constraints
                initial_parm_vals = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., 0.],
                ndim = 12,
                nwalkers = 120, # ndim * 10
                niterations = 1000, # number of steps in mcmc sampler
                interval = 10, # the interval over which we save the data.
                out_tag = '',
                out_dir = './'
              ):
    
    # initialize the likelihood function
    logposterior = UHE_likelihood(UHE_fluence, data)
    print 'Initial Likelihood'
    logposterior(initial_parm_vals)

    # Set up the emcee Ensemble sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, logposterior, threads=1)
    
    # Dither the initial positions
    pos = np.array( [initial_parm_vals*(1. + 1.e-2*np.random.randn(ndim)) for i in range(nwalkers)] )
    # Add a little bit of composition
    for k in range(0,len(pos)):
        pos[k][3:7] += np.random.uniform(0.,0.01,4)

    # Set up an interval of iterations over which to output the chains.
    N_loops = niterations / interval
    print interval, N_loops, N_loops*interval, niterations
    for N in range(0,N_loops):
        print 'N=', N
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
    

