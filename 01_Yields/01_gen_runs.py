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


#log_energy_range   = np.arange(17.,25.1,0.1) # in log eV
#log_redshift_range = np.arange(-4., 1.1, 0.1)
#Z_list = [1,2,7,14, 26] 
#A_list = [1,4,14,28,56] 

log_energy_range   = np.arange(inputs['log10_energy_min'],   inputs['log10_energy_max']+inputs['log10_energy_step'],    inputs['log10_energy_step']) # in log eV
log_redshift_range = np.arange(inputs['log10_redshift_min'], inputs['log10_redshift_max']+inputs['log10_redshift_step'],inputs['log10_redshift_step'])
Z_list = inputs['Z_list']
A_list = inputs['A_list']

print 'Number of CRPropa runs:', len(Z_list)*len(log_energy_range)*len(log_redshift_range)

#outdir_base = '/halo_nobackup/lenssim/romerowo/crpropa_production_run_outputs/20160506/'
#inputdir_base = './inputs'

outdir_base   = inputs['CRPropa_output_dir']
inputdir_base = inputs['CRPropa_inputs_dir']
crane_dir     = os.environ['CRANEDIR']

global_count=0
for k in range(0,len(Z_list)):
      print 'Setting up runs for Z=',Z_list[k]
      inputdir = inputdir_base + '/Z_%d/'%Z_list[k]
      outdir = outdir_base + '/Z_%d/'%Z_list[k]
      if(not os.path.exists(inputdir)):
            print 'Creating input directory', inputdir
            os.makedirs(inputdir)
      if(not os.path.exists(outdir)):
            print 'Creating output directory', outdir
            os.makedirs(outdir)
      count = 0
      for e in log_energy_range:
          for z in log_redshift_range:
              if np.abs(z)<1.e-6:
                  z = 0.
              #print Z_list[k], A_list[k], e, z
              f = open(inputdir+'input.%d'%count,'w')
              com = '%s/crpropa_run.py -E %1.1f -Z %d -A %d -log10_redshift %1.1f -od %s -crane_dir %s'%(inputs['CRANE_dir']+'/01_Yields',e, Z_list[k], A_list[k], z, outdir, os.environ['CRANEDIR'])
              #print count, com
              f.write(com)
              f.close()
		      #print com
              count += 1
              global_count += 1
      # edit pbs_file for computing cluser run
      tmp_file = open('tmp.sh', 'w')
      # print 'count', count
      for line in file('pbs_run_Z_%d.sh'%Z_list[k]):
          #print line[:-1]
          new_line = line
          if ('#PBS -J' in line and ('\'PBS -J\'' not in line and '\'#PBS -J\'' not in line and '\"#PBS -J\"' not in line)):
            new_line = '#PBS -J 0-%d\n'%(count-1)
          if('export PBS_O_WORKDIR' in line and line[0]!='#'):
            new_line = 'export PBS_O_WORKDIR=%s/Z_%d\n'%(outdir_base,Z_list[k])
          if('export INPUT_DIR' in line and line[0]!='#'):
            new_line = 'export INPUT_DIR=%s/Z_%d\n'%(inputdir_base,Z_list[k])
          if('$PBS_ARRAY_INDEX' in line and line[0]!='#' and 'mv' not in line):
            new_line = '%s/utils/run_python_command.py $INPUT_DIR/input.$PBS_ARRAY_INDEX > output.$PBS_ARRAY_INDEX\n'%(crane_dir)

          tmp_file.write(new_line)
      tmp_file.close()
      os.system('mv tmp.sh pbs_run_Z_%d.sh'%Z_list[k])
#print 'global_count', global_count

