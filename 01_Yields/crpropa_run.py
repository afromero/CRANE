#!/usr/bin/env python
'''
2016 May 6
ARW, JPL/Caltech
CRpropa based UHECR and neutrino propagation simulations
'''
import sys
if __name__ == "__main__":
      import argparse
      parser=argparse.ArgumentParser(description='CRpropa-based simulation for nuclei and neutrinos')

      parser.add_argument("-E",              "--log10_energy",       default=20., 	help="Nuclear Energy in log10 eV basis", type=float)
      parser.add_argument("-Z",              "--charge",       default=1, 	      help="Nuclear Charge", type=int)
      parser.add_argument("-A",              "--mass",         default=1, 	      help="Nuclear Mass", type=int)
      parser.add_argument("-log10_redshift", "--log10_redshift",     default=0.001826, help="redshift in log10 basis", type=float)
      parser.add_argument("-nparticles",     "--num_particles",default=10000, 	help="number of particles to simulate", type=int)
      parser.add_argument("-o",              "--output_tag",   default=None, 	help="Output tag, in quotes",type=str)
      parser.add_argument("-od",             "--outputdir",    default='./', 	help="Output directory, in quotes",type=str)
      parser.add_argument("-crane_dir",      "--crane_dir",    default='../', 	help="Output directory, in quotes",type=str)
      #KLUDGE TO GET ARGPARSE TO READ NEGATIVE VALUES
      for i, arg in enumerate(sys.argv):
            if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg
      #PAESE ARGUMENTS
      args=parser.parse_args()
      # ensure the output director ends with '/' to keep the format consisent.
      if(args.outputdir[-1]!='/'):
	      args.outputdir = args.outputdir+'/'
      if(args.output_tag==None):
            args.output_tag = 'run_E_%1.1f_Z_%d_A_%d_r_%1.1f'%(args.log10_energy, args.charge, args.mass, args.log10_redshift)
            print args.output_tag

      print ''
      print 'SIMULATION PARAMETERS' 
      print args
      print ''


      import os
      inputs = eval(open("%s/inputs.txt"%args.crane_dir).read()) 
      #print 'inputs'
      #print inputs

      # point to the CRPropa installation dir
      sys.path.insert(0, inputs['CRPropa_dir'])
      from crpropa import *

      neutrinos = True
      photons = False

      # initialize CRpropa module list
      m = ModuleList()

      # the arguments are the minimum and maximum step size.
      m.add(SimplePropagation())
      #m.add(SimplePropagation(10*kpc, 10*Mpc))

      #  update redshift and apply adabiatic loss according to travelled distance
      m.add(Redshift())

      # add interactions due to CMB, IRB
      # the full list is
      # CMB, IRB, IRB_Kneiske04, IRB_Stecker05, IRB_Franceschini08, IRB_Finke10, IRB_Dominguez11, 
      # IRB_Gilmore12, IRB_withRedshift_Kneiske04, IRB_withRedshift_Stecker05, 
      # IRB_withRedshift_Franceschini08, IRB_withRedshift_Finke10, IRB_withRedshift_Dominguez11, 
      # IRB_withRedshift_Gilmore12, URB_Protheroe96 
      m.add(PhotoPionProduction(CMB, photons, neutrinos))
      m.add(PhotoPionProduction(IRB, photons, neutrinos))
      m.add(PhotoDisintegration(CMB)) 
      m.add(PhotoDisintegration(IRB))
      m.add(NuclearDecay(photons, neutrinos))
      m.add(ElectronPairProduction(CMB))
      m.add(ElectronPairProduction(IRB))

      m.add(MinimumEnergy(10**15 * eV))

      # observer for hadrons
      obs1 = Observer()
      obs1.add(ObserverPoint())
      obs1.add(ObserverNeutrinoVeto())  # we don't want neutrinos here
      nucleon_out_filename = args.outputdir+'nucleons_'+args.output_tag+'.txt' 
      #print nucleon_out_filename
      output1 = TextOutput(nucleon_out_filename, Output.Event1D)
      obs1.onDetection( output1 )
      m.add(obs1)
      # observer for neutrinos
      obs2 = Observer()
      obs2.add(ObserverPoint())
      obs2.add(ObserverNucleusVeto())  # we don't want hadrons here
      neutrino_out_filename = args.outputdir+'neutrinos_'+args.output_tag+'.txt' 
      #print neutrino_out_filename
      #exit()
      output2 = TextOutput(neutrino_out_filename, Output.Event1D)
      obs2.onDetection( output2 )
      m.add(obs2)

      import time
      start = time.clock()
      # source: protons with power-law spectrum from uniformly distributed sources with redshift z = 0-3
      source = Source()

      # source at a given location
      #source.add( SourcePosition(200 * Mpc) )
      # source at a given redshift
      source.add( SourcePosition( redshift2ComovingDistance ( 10**(args.log10_redshift) ) ) )

      # set the source redshift according to the distance to 0.
      source.add(SourceRedshift1D())

      # Set source energy
      source.add( SourceEnergy(10**(args.log10_energy) * eV) )
      #source.add(SourcePowerLawSpectrum(10**18 * eV, 10**25 * eV, -1))

      # set nuclear species
      source.add( SourceParticleType(nucleusId(args.mass, args.charge)) )
      # source.add(SourceParticleType(nucleusId(1, 1)))

      print source.getDescription()
      candidate = source.getCandidate()
      print candidate.getRedshift()

      # run simulation for 1000 primaries and propagate all secondaries
      m.setShowProgress(True)
      #num_particles = 10000
      m.run(source, args.num_particles, True)
      print 'time elapsed',  time.clock() - start
      output1.close()
      output2.close()

      # REDUCE TO HISTOGRAMS
      import numpy as np

      # get hadron data
      data = np.genfromtxt(nucleon_out_filename, names=True)




      # set energy array for plotting
      logE0 = np.log10(data['E0']) + 18
      logE  = np.log10(data['E']) + 18
      print 'Number of events', len(data)
      # get arrival nuclear charges and atomic mass numbers
      count=0
      Z = np.array([chargeNumber(id) for id in data['ID'].astype(int)])  # element
      A = np.array([massNumber(id) for id in data['ID'].astype(int)])  # atomic mass number

      # index data by atomic mass number
      u, counts = np.unique(A, return_counts=True)
      #print u

      nbins = 100
      rng_min = 15.
      rng_max = 25.
      h_list = []
      for a in u:
          h, e = np.histogram(logE[A==a], bins=nbins, range=(rng_min, rng_max))
          h_list.append(h)
          print a, np.sum(h)
      print e

      for k in range(0,len(u)):
            print u[k], np.sum(h_list[k])

      if(os.path.getsize(neutrino_out_filename)!=0):
            d = np.genfromtxt(neutrino_out_filename, names=True)
            # histogram neutrino counts by energy.
            # note, the data does tabulate neutrino flavors.
            h_nu, e = np.histogram(np.log10(d['E']) + 18, bins=nbins, range=(rng_min, rng_max))
            print np.sum(h),h,e
            #hist(log10(d['E']) + 18, bins=nbins, range=(rng_min, rng_max), histtype='step', normed=False, 
            #      log=True, label='At z=0, Neutrinos')
      if(os.path.getsize(neutrino_out_filename)==0):
            print 'NO NEUTRINOS PRODUCED!'
            h_nu = 0*h_list[0]

      np.savez(args.outputdir+args.output_tag+'.npz', input_energy =  args.log10_energy,
                                                      num_particles = args.num_particles, 
                                                      A_list = u, 
                                                      energy_bins = e, 
                                                      nucleon_histograms=h_list,
                                                      neutrino_histogram = h_nu)
      os.system('rm %s'%(nucleon_out_filename))
      os.system('rm %s'%(neutrino_out_filename))


      exit()

      dat = np.load(args.outputdir+args.output_tag+'.npz')
      for k, v in dat.items():
            print '\t',k,v

      import matplotlib
      matplotlib.use('Agg')
      from pylab import *
      import numpy as np
      rcParams['font.size'] = 16

      plt.figure(figsize=(10, 7))
      ax = subplot(111)
      ax.set_yscale('log')
      # histogram input particle hadrons
      print 'dat[\'A_list\']', dat['A_list'] 
      print 'dat[\'A_list\'].size',  dat['A_list'].size
      print 'len(dat[\'A_list\'])', len(dat['A_list'])
      bin_width = dat['energy_bins'][1] - dat['energy_bins'][0]
      for k in range(0, dat['A_list'].size):
            print k
            #bar(dat['energy_bins'][:-1], dat['nucleon_histograms'][k], (rng_max-rng_min)/float(nbins),  label='A=%d'%dat['A_list'][k], alpha=0.5)
            step(dat['energy_bins'][:-1] + bin_width , dat['nucleon_histograms'][k],  label='A=%d'%dat['A_list'][k], alpha=1)

      step(dat['energy_bins'][:-1]+bin_width, dat['neutrino_histogram'],  color='k', label='neutrinos', alpha=1., lw=2)
      xlim(15, 25.)
      xlabel(r'$\log_{10}(E/{\rm eV})$')
      ylabel(r'Number of Particles')
      plt.title('%d particles of , Z=%d, A=%d, log10(Energy / eV) = %1.1f , redshift=%1.4f'%(args.num_particles, args.charge, args.mass, args.log10_energy, 10**(args.log10_redshift)), fontsize=14)
      plt.legend(loc = 0, fontsize=10)
      plt.xticks([15.,16.,17.,18.,19.,20.,21., 22., 23., 24., 25.])
      grid(True, which='both')
      plt.savefig(args.outputdir + args.output_tag+'.png')

      dat.close()
      exit()

      ############################################################################################################
      ############################################################################################################
      ############################################################################################################

