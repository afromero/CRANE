import corner
import pickle
import os
from pylab import *

inputs = eval(open("../inputs.txt").read())
sys.path.insert(0, inputs['CRANE_dir']+'/02_Inference')


def load_chains_and_probs(tag='', dirc='.'):
   
  for k in range(0,1000):
    chain_fname  = dirc+'/uhe_emcee_chain_'  + tag +'_%d.p'%k
    lnprob_fname = dirc+'/uhe_emcee_lnprob_' + tag +'_%d.p'%k
    #print chain_fname
    if(os.path.os.path.isfile(chain_fname)):
        if(k==0):
            big_chain = pickle.load( open( chain_fname, "rb" ) )
            #print big_chain.shape
            big_lnprob = pickle.load( open( lnprob_fname, "rb" ) )
            #print big_lnprob.shape
        if(k>0):
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
chain, lnprob = load_chains_and_probs(tag = 'run01', dirc = '/halo_nobackup/eva/romerowo/crane_inference_outputs/20161103')

nwalkers, niterations, ndim = chain.shape

n_start = int(0.5*float(niterations))
samples = chain[:, n_start:, :].reshape((-1, ndim))

print 'PLOTTING FIGURES'
print '\tplotting likelihood'
figure()
subplot(211)
for k in range(0,nwalkers):
    plot(lnprob[k,n_start:], 'k', alpha=0.2)
subplot(212)
hist(np.ravel(lnprob[:,n_start:]))

#print '\tplotting E_max likelihood'
#figure()
#plot(chain[:,:,2], lnprob[:,:], 'k.')

print '\tplotting fluences'
figure(figsize=(8,12))
mx_nu = 1.e-308
mx_cr = 1.e-308
colors = cm.rainbow(np.linspace(0, 1, 100))
for k in range(0,1000):
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



print '\tplotting corner'
labels = ['norm','spec_ind', 'E_max', 'f_He', 'f_N', 'f_Si', 'f_Fe', 'src_ind', 'z1', 'z2', 'z3', 'ln(beta)']
print 'samples.shape', samples.shape
samples[:,11]= np.log(samples[:,11])
fig = corner.corner(samples, labels=labels)
savefig('marginalized_paramters_distributions.png')
#fig.set_size_inches(8, 8)
show()
