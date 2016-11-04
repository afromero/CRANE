from pylab import *
import os

# Script to check that the runs completed, find runs that did not, and try to fix them
if('CRANEDIR' not in os.environ.keys()):
    print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print 'The \'CRANEDIR\' environemnt variable has not been set'
    print '\tFor example, in csh'
    print '\tsetenv CRANEDIR /home/romerowo/CRANE'
    print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    exit()
# load inputs
inputs = eval(open("%s/inputs.txt"%os.environ['CRANEDIR']).read()) 

# LOOP THROUGH INPUTS
def parse_input_file(fnm):
    run_index = fnm.split('.')[1]
    for line in file(fnm):
        #print line
        out_dir = line.split('-od')[1].split()[0]
        out_name = 'run_E_%s_Z_%s_A_%s_r_%s'%(line.split('-E')[1].split()[0], line.split('-Z')[1].split()[0], line.split('-A')[1].split()[0], line.split('-log10_redshift')[1].split()[0])
        out_fnm = out_dir + out_name + '.npz'
    return run_index, out_fnm

def check_finished(output_fnm):
    for line in file(output_fnm):
        if('[ Finished ]') in line:
            return True
    return False

def create_pbs(run_index):
    new_file = open('pbs_run_Z_26_%s.sh'%run_index,'w')
    for line in file('pbs_run_Z_26_X.sh'):
       new_line = line
       if('.X' in line):
         new_line = line.replace('.X', '.%s'%run_index)
         #print new_line
       new_file.write(new_line)
    new_file.close()
    #os.system('qsub pbs_run_Z_26_%s.sh'%run_index)
counter = 0
for Z in inputs['Z_list']:
  if(Z!=26): continue
  dirc = inputs['CRPropa_inputs_dir'] + '/Z_%d'%Z
  files = os.listdir(dirc)
  for fnm in files:
    #print fnm
    run_index, out_fnm = parse_input_file(dirc + '/' + fnm)
    #print out_fnm
    if( not os.path.isfile(out_fnm)):
       print out_fnm
       counter += 1
       # check if CRPropa text output files are present
       nu_fnm = out_fnm.replace('.npz', '.txt').replace('/run_', '/neutrinos_run_')
       nuc_fnm = out_fnm.replace('.npz', '.txt').replace('/run_', '/nucleons_run_')
       output_fnm = inputs['CRPropa_output_dir'] + '/Z_%d'%Z + '/output.%s'%run_index
       print '\t', run_index
       if not os.path.isfile(output_fnm):
            print 'No output file found!', output_fnm
            continue
       print '\t', 'Finished:', check_finished(output_fnm)
       if(not os.path.isfile(nuc_fnm)): 
           print '\tNucleon File exists'
           exit()
       if(not os.path.isfile(nu_fnm)): 
           print '\tNeutrino File exists'
           exit()
       create_pbs(run_index)
       #print nu_fnm, nuc_fnm
print 'Files Missing:', counter

    # check that the output file exists and its size is more than zero
    

