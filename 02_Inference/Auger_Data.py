import numpy as np

# Just for reference from PRL, errors are gaussian
d_log10_E=0.1
log10_E=[ 18.65, 18.75, 18.85, 18.95, 19.05, 19.15, 19.25, 19.35, 19.45, 19.55, 19.65, 19.75, 19.85, 19.95, 20.05 ] # these are center bin values
nEvents = [ 4778, 3159, 2162, 1483, 1052, 699, 451, 323, 200, 110, 43 ,28 , 23 ,5 , 2 ]
#Exposure_Auger = 12790.*np.ones(len(log10_E)) # km2 sr yr
Exposure = 12790. * 3.1556926e17 * np.ones(len(log10_E)) # cm2 sr s

# we should add uncertainties in the Energy and Exposure estimate.

# Need to add mean Xmax and rms Xmax to this file

