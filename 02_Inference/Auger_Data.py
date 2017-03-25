import numpy as np
from scipy.special import factorial

d_log10_E=0.1
log10_Energy_low_edges = np.arange(18.6, 25. ,0.1)
Counts = np.zeros(len(log10_Energy_low_edges))
Exposure = np.zeros(len(log10_Energy_low_edges))

# The values below are from
#"Measurement of the cosmic ray energy spectrum using hybrid events of the Pierre Auger Observatory," M. Settimo for the Pierre Auger Collaboration, Eur. Phys. J. Plus 127 (2012) 87 [doi:10.1140/epjp/i2012-12087-9] [arXiv:1208.6574]
# these are center bin energy values
#log10_E = np.array([ 18.65, 18.75, 18.85, 18.95, 19.05, 19.15, 19.25, 19.35, 19.45, 19.55, 19.65, 19.75, 19.85, 19.95, 20.05 ]) 
# nEvents = np.array([ 4778, 3159, 2162, 1483, 1052, 699, 451, 323, 200, 110, 43 ,28 , 23 ,5 , 2 ])
#Constant_Exposure = 12790. * 3.1556926e17 # cm2 sr s 
#Exposure =  Constant_Exposure * np.ones(len(log10_Energy_low_edges)) 

# The values below are from:
# "The flux of ultra-high energy cosmic rays after ten years of operation of the Pierre Auger Observatory," 
# I. Valino for the Pierre Auger Collaboration, The Pierre Auger Observatory: Contributions to the 34thInternational Cosmic Ray Conference (ICRC 2015)[arXiv:1509.03732]
# these are center bin energy values
#log10_E = np.array([ 18.65, 18.75,   18.85,   18.95,   19.05,   19.15,   19.25,   19.35,   19.45,   19.55,   19.65,   19.75,   19.85,   19.95,   20.05,   20.15 ]) 
nEvents = np.array([  26325, 16317 ,  10991,   7585,    5221,    3491,    2238,    1405,    888,     569,     267,     130,     54,      10,      4,       1])
# expaVals extracted from digitizing Figure 3 of the reference. Values are in km^2 yr sr
expVals = np.array([60557.5, 57222.9, 55609.2, 55366.3, 55896.3, 52759.4, 55063.0, 54467.3, 55748.4, 53567.6, 56308.1, 58779.1, 58126.4, 58488.4, 43708.9, 43147.9])

Counts[0:len(nEvents)]   = nEvents
Exposure[0:len(nEvents)] = expVals
Exposure[len(nEvents)-1:] = Exposure[len(nEvents)-1] # set the exposure at higher energies to the last non-zero value
Exposure *= 3.1556926e17 # cm2 sr s 

# we should add uncertainties in the Energy and Exposure estimate. These are not published.

# Need to add mean Xmax and rms Xmax to this file
X_max_log10_E_low_edge   = np.arange(17.8, 19.51, d_log10_E)
X_max_log10_E_upper_edge = np.arange(17.9, 19.61, d_log10_E)
X_max_Num       = np.array([3667,  3365,  2859,  2436,  1984,  1442,  1150,  832,   591,   431,   324,   246,   174,   129,   96,    64,    44,    38]) 
X_max_Mean      = np.array([714.3, 723.3, 730.6, 740.5, 745.4, 752.0, 757.2, 757.3, 759.2, 758.1, 761.7, 769.1, 763.7, 768.6, 775.9, 777.3, 783.2, 774.1]) 
X_max_Mean_err  = np.array([1.5,   1.7,   1.6,   2.0,   1.9,   2.3,   2.5,   2.4,   2.9,   2.6,   2.5,   3.9,   3.1,   4.0,   4.8,   6.2,   9.2,   4.7]) 
X_max_Mean_sysU = np.array([9.88,  9.93,  10.0,  10.1,  10.2,  10.3,  10.4,  10.6,  10.7,  10.9,  11.1,  11.3,  11.5,  11.7,  11.9,  12.1,  12.3,  12.6])
X_max_Mean_sysD = np.array([7.8,   7.8,   7.8,   7.8,   7.8,   7.9,   8.0,   8.1,   8.2,   8.3,   8.5,   8.6,   8.8,   8.9,   9.1,   9.2,   9.3,   9.6])
X_max_RMS       = np.array([55.8,  60.0, 62.9,   63.4,  66.6,  62.4,  58.6,  56.3,  57.9,  45.2,  41.7,  50.1,  38.0,  44.2,  40.4,  46.9,  48.7,  24.4])
X_max_RMS_err   = np.array([2.8,   3.1,  2.9,    3.4,   3.0,   3.9,   4.4,   3.6,   4.5,   3.3,   3.0,   3.5,   4.0,   4.9,   6.3,   7.0,   12.0,  5.5])
X_max_RMS_sysU  = np.array([5.9,   5.5,  5.3,    5.2,   5.1,   5.1,   5.1,   5.1,   5.1,   5.1,   5.1,   5.1,   5.1,   5.1,   5.1,   5.0,   5.0,   5.0])
X_max_RMS_sysD  = np.array([6.4,   5.7,  5.4,    5.3,   5.2,   5.1,   5.1,   5.1,   5.1,   5.1,   5.1,   5.1,   5.1,   5.1,   5.1,   5.1,   5.1,   5.4])
'''
 17.8 17.9 17.850  3667 714.3 1.5 9.88 7.8 55.8 2.8 5.9 6.4
 17.9 18.0 17.949  3365 723.3 1.7 9.93 7.8 60.0 3.1 5.5 5.7
 18.0 18.1 18.048  2859 730.6 1.6 10.0 7.8 62.9 2.9 5.3 5.4
 18.1 18.2 18.148  2436 740.5 2.0 10.1 7.8 63.4 3.4 5.2 5.3
 18.2 18.3 18.247  1984 745.4 1.9 10.2 7.8 66.6 3.0 5.1 5.2
 18.3 18.4 18.348  1442 752.0 2.3 10.3 7.9 62.4 3.9 5.1 5.1
 18.4 18.5 18.448  1150 757.2 2.5 10.4 8.0 58.6 4.4 5.1 5.1
 18.5 18.6 18.549   832 757.3 2.4 10.6 8.1 56.3 3.6 5.1 5.1
 18.6 18.7 18.647   591 759.2 2.9 10.7 8.2 57.9 4.5 5.1 5.1
 18.7 18.8 18.747   431 758.1 2.6 10.9 8.3 45.2 3.3 5.1 5.1
 18.8 18.9 18.849   324 761.7 2.5 11.1 8.5 41.7 3.0 5.1 5.1
 18.9 19.0 18.947   246 769.1 3.0 11.3 8.6 50.1 3.5 5.1 5.1
 19.0 19.1 19.048   174 763.7 3.1 11.5 8.8 38.0 4.0 5.1 5.1
 19.1 19.2 19.146   129 768.6 4.0 11.7 8.9 44.2 4.9 5.1 5.1
 19.2 19.3 19.245    96 775.9 4.8 11.9 9.1 40.4 6.3 5.1 5.1
 19.3 19.4 19.340    64 777.3 6.2 12.1 9.2 46.9 7.0 5.0 5.1
 19.4 19.5 19.446    44 783.2 9.2 12.3 9.3 48.7 12. 5.0 5.1
 19.5 20.0 19.620    38 774.1 4.7 12.6 9.6 24.4 5.5 5.2 5.4
'''



# This function models ther number of events given an input fluence curve
def Counts_Model(nuclear_fluence, log10_E_shift = 0.):
    #print 'len check', len(Exposure), len(nuclear_fluence)
    arrival_Counts = nuclear_fluence * Exposure[-len(nuclear_fluence):] * np.log(10.) * d_log10_E
    if(log10_E_shift == 0.):
        return arrival_Counts

    res       = log10_E_shift%d_log10_E
    bin_shift = int(round((log10_E_shift - res)/d_log10_E))
    bin_res = res/d_log10_E
    # shift and re-weight the bins
    observed_Counts = np.roll(arrival_Counts,bin_shift)
    observed_Counts = (1.-bin_res)*np.roll(arrival_Counts,bin_shift) + (bin_res)*np.roll(arrival_Counts,bin_shift+1)
    if(bin_shift>=0):
        observed_Counts[0:bin_shift+1] *= 0
    if(bin_shift<0):
        observed_Counts[bin_shift:] *= 0
    return observed_Counts


# This is useful to have as a pre-defined array for the Poisson likelihood
log_factorial = []
for k in range(0,len(Counts)):
    if np.isfinite(factorial(Counts[k])):
        log_factorial.append(np.log(factorial(Counts[k])))
    else:
        # use Stirling's approximation 
        # (Note: there is an approximation by Ramanujan that is more accurate, should implement in the future.)
        val = Counts[k]*np.log(Counts[k]) - Counts[k]
        log_factorial.append(val)
log_factorial = np.array(log_factorial)

