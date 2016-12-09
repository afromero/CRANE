import corner
import pickle
import os
from pylab import *

inputs = eval(open("../inputs.txt").read())
sys.path.insert(0, inputs['CRANE_dir']+'/02_Inference')


def load_chains_and_probs(tag='', dirc='.', k_start=0, k_end=1000):
  for k in range(k_start,k_end):
    chain_fname  = dirc+'/uhe_emcee_chain_'  + tag +'_%d.p'%k
    lnprob_fname = dirc+'/uhe_emcee_lnprob_' + tag +'_%d.p'%k
    #print chain_fname
    if(os.path.os.path.isfile(chain_fname)):
        if(k==k_start):
            big_chain = pickle.load( open( chain_fname, "rb" ) )
            #print big_chain.shape
            big_lnprob = pickle.load( open( lnprob_fname, "rb" ) )
            #print big_lnprob.shape
        if(k>k_start):
            chain = pickle.load( open( chain_fname, "rb" ) )
            #print chain.shape
            lnprob = pickle.load( open( lnprob_fname, "rb" ) )
            #print lnprob.shape
            big_chain  = np.concatenate([big_chain,  chain],  axis=1)
            big_lnprob = np.concatenate([big_lnprob, lnprob], axis=1)
        #print big_chain.shape
        #print big_lnprob.shape
  return big_chain, big_lnprob

print 'LOADING FLUENCE CALCULATOR'
from Fluence import *
# read inputs

# load npz file (NOTE: this is a large file)
print 'Loading Lookup Table'
f = np.load(inputs['Lookup_Table'])
print 'Finished Loading'
input_Z = f['input_Z']
input_A = f['input_A']
input_log10_E = f['input_log10_E']
input_log10_z = f['input_log10_z']

observed_log10_E = f['observed_log10_E']
observed_A = f['observed_A']
uhe_LUT = f['LUT']

fcalc = UHE_fluence(input_Z, input_A, input_log10_E, input_log10_z, observed_A, observed_log10_E, uhe_LUT)
        
print 'LOADING CHAINS'
#this one used the Auger data published in 2012. No X_max was included in the likelihood. Only counts
#chain, lnprob = load_chains_and_probs(tag = 'run01', dirc = '/halo_nobackup/eva/romerowo/crane_inference_outputs/20161103', k_start=500, k_end=668)

#this one used the Auger data published in 2015 (one decade). X_max was also included in the likelihood. MCMC was initialized with various points around the liklihood space.
#chain, lnprob = load_chains_and_probs(tag = 'run02', dirc = '/halo_nobackup/eva/romerowo/crane_inference_outputs/20161110', k_start=0, k_end=200)

#this one used the Auger data published in 2015 (one decade). X_max was also included in the likelihood. MCMC was initialized around what appears to be a well-defined global minimum.
# chain, lnprob = load_chains_and_probs(tag = 'run03', dirc = '/halo_nobackup/eva/romerowo/crane_inference_outputs/20161114', k_start=0, k_end=200)

#this one included priors on astrophysical distribution of sources.
chain, lnprob = load_chains_and_probs(tag = 'run04', dirc = '/halo_nobackup/eva/romerowo/crane_inference_outputs/20161128', k_start=0, k_end=200)

#print lnprob.shape
#max_LL =  np.max(lnprob)
#print np.where(lnprob>0.99*max_LL)
#exit()

nwalkers, niterations, ndim = chain.shape

n_start = int(0.5*float(niterations))
samples = chain[:, n_start:, :].reshape((-1, ndim))


print 'PLOTTING FIGURES'
print '\tplotting likelihood'
figure()
subplot(211)
for k in range(0,nwalkers):
    plot(lnprob[k,n_start:], 'k', alpha=0.2)
ylim(-300.,0.)
subplot(212)
hist(np.ravel(lnprob[:,n_start:])[np.isfinite(np.ravel(lnprob[:,n_start:]))], bins=np.arange(-300.,0., 1.))
#show()
#print '\tplotting E_max likelihood'
#figure()
#plot(chain[:,:,2], lnprob[:,:], 'k.')

print '\tplotting fluences'
figure(figsize=(8,12))
mx_nu = 1.e-308
mx_cr = 1.e-308
colors = cm.rainbow(np.linspace(0, 1, 100))
for k in range(0,100):
    if(k%10==0): print '\t',k
    i = np.random.randint(0,len(samples))
    parms = samples[i,:]
    fluences = fcalc.fluence_model(*parms)
    ax=subplot(211)
    ax.set_yscale('log')
    plot(observed_log10_E, fluences[0,:], 'b-', alpha=0.1)
    ylabel('Fluence, cm$^{-2}$ s$^{-1}$ sr$^{-1}$')
    xlabel('log$_{10}$(Energy / eV)')
    grid(True)
    title('Neutrino Fluence')
    mx_nu = np.max([mx_nu, np.max(fluences[0,:])])
    #y1,y2 = ax.get_ylim()
    ylim(1.e-10*mx_nu, mx_nu)
    xlim(15., 25.)
    #legend(loc=0)

    ax2=subplot(212)
    ax2.set_yscale('log')
    #nuc_fluence = np.sum(fluences[1:], axis=0)
    ax2.plot(observed_log10_E, np.sum(fluences[1:], axis=0), 'b-', lw=1, alpha=0.1,  label='A=1-56')
    #legend(loc=3)
    xlim(15., 25.)
    mx_cr = np.max([mx_cr, np.max(np.sum(fluences[1:], axis=0))])
    ylim(1.e-20*mx_cr, mx_cr)
    y1, y2 = ax2.get_ylim()
    fill_between([15., 18.9], [y1, y1], [y2,y2], facecolor='none', hatch='//', edgecolor='gray', linewidth=0.0)
    #legend(loc=0)
    ylabel('Fluence, cm$^{-2}$ s$^{-1}$ sr$^{-1}$')
    title('Nuclear Fluence')
    xlabel('log$_{10}$(Energy / eV)')
    subplots_adjust(hspace=0.4)
    grid(True)

savefig('inferred_fluences.png')
##########
#print 'showing'
#show()

#levels = 1.0 - np.exp(-0.5 * np.arange(0.5, 3.1, 0.5) ** 2)
levels = 1.0 - np.exp(-0.5 * np.arange(1., 3.1, 1.) ** 2)
print 'levels', levels


print '\tplotting corner'
labels = ['norm','spec_ind', 'E_max', 'f_He', 'f_N', 'f_Si', 'f_Fe', 'src_ind', 'z1', 'z2', 'z3', 'ln(beta)']
print 'samples.shape', samples.shape
samples[:,11]= np.log(samples[:,11])
hist2d_kwargs = {'plot_datapoints' : False, 'levels' : levels}
fig = corner.corner(samples, 
                    labels=labels, 
                    smooth = 1., 
                    smooth1d = 1., 
                    color='k', fill_contours=True, 
                    **hist2d_kwargs)
savefig('marginalized_paramters_distributions.png')
#fig.set_size_inches(8, 8)
show()
