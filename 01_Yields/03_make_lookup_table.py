from pylab import *
import os

if('CRANEDIR' not in os.environ.keys()):
    print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print 'The \'CRANEDIR\' environemnt variable has not been set'
    print '\tFor example, in csh'
    print '\tsetenv CRANEDIR /home/romerowo/CRANE'
    print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    exit()

inputs = eval(open("%s/inputs.txt"%os.environ['CRANEDIR']).read()) 

log_energy_range   = np.arange(inputs['log10_energy_min'],   inputs['log10_energy_max']+inputs['log10_energy_step'],    inputs['log10_energy_step']) # in log eV
log_redshift_range = np.arange(inputs['log10_redshift_min'], inputs['log10_redshift_max']+inputs['log10_redshift_step'],inputs['log10_redshift_step'])
Z_list = inputs['Z_list']
A_list = inputs['A_list']

A_observed          = range(0,57) # atomic mass at Earth, 0 is neutrinos, 1 and higher 
log_energy_observed = np.arange(15.,25.0,0.1) # energy at Earth

print 'initializing large LUT'
LUT = np.zeros((len(Z_list), len(log_energy_range), len(log_redshift_range), len(A_observed), len(log_energy_observed)))
print 'LUT.shape',LUT.shape

# Start loop to read through the outputs

for Z_index in range(0,len(Z_list)):
  print 'Nuclear Species (Z,A):',(Z_list[Z_index], A_list[Z_index])
  for e_index in range(0,len(log_energy_range)):
    if(e_index%10==0): print '\tSource Energy:', log_energy_range[e_index]
    for log_10_z_index in range(0,len(log_redshift_range)):
      # get file name of the CRPropa outputs
      fnm = inputs['CRPropa_output_dir'] + '/Z_%d/run_E_%1.1f_Z_%d_A_%d_r_%1.1f.npz'%(Z_list[Z_index],log_energy_range[e_index],Z_list[Z_index],A_list[Z_index],log_redshift_range[log_10_z_index])
      if(not os.path.exists(fnm)):
        print '\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        print 'WARNING:'
        print 'File not found: %s'%fnm
        print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n'
        continue
      # load the npz file with the outputs
      dat = np.load(fnm)

      # fill neutrino spectra in the entries corresponding to a_index = 0
      LUT[Z_index, e_index, log_10_z_index, 0, : ] = dat['neutrino_histogram']/1.e4
    
      # loop through the output nuclear masses and fill in the LUT
      a_index=0
      for a_out in dat['A_list']:
        if(a_out != dat['A_list'][a_index]):
            print 'a_out, a_index, dat[\'A_list\'][a_index]', a_out, a_index, dat['A_list'][a_index]
            exit()
        LUT[Z_index, e_index, log_10_z_index, a_out, : ] = dat['nucleon_histograms'][a_index]/1.e4
        a_index+=1

      # close npz file with the CRPRopa outputs
      dat.close()

np.savez(inputs['Lookup_Table'], input_Z          = Z_list, 
                                 input_A          = A_list,
                                 input_log10_E    = log_energy_range,
                                 input_log10_z    = log_redshift_range,
                                 observed_A       = A_observed,
                                 observed_log10_E = log_energy_observed,
                                 LUT =  LUT)


exit()
