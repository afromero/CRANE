import numpy as np
from scipy.interpolate import Akima1DInterpolator

class Xmax_modeler:
    def __init__(self, log10_Energy_array):
        print 'Initializing Xmax_modeler'
        self.log10_Energy_array = log10_Energy_array
        # initiate modeling parameters derived from Conex
        # model names
        self.model = ["qgsjetII04" ,"eposlhc"] 
        # The nuclear masses used to obtain Conex parameters
        self.A_Conex = [1,4,14,28,40,56] # p, He, N, Si, Ca, Fe
        self.A_array = range(1,57)

        # Define 2x6 arrays. The first index is for the model, and the second is for A_Conex.
        self.Offset_A_Conex = np.array([[ 652.6, 660.2, 652.8, 650.7, 646.0, 642.5],
                                        [ 674.3, 687.8, 685.1, 678.2, 673.7, 669.1]])
        self.Offset_B_Conex = np.array([[ 63.1,  60.8,  59.6,  59.3,  59.3,  59.6 ],
                                        [ 63.7,  64.1,  62.8,  63.1,  63.1,  63.3 ]])
        self.Tau_A_Conex    = np.array([[ 27.5,  19.4,  14.7,  12.1,  11.2,  10.3 ],
                                        [ 26.1,  17.2,  11.5,   9.5,   8.6,   7.8 ]])
        self.Tau_B_Conex    = np.array([[ -1.7,  -1.2,  -0.8,  -0.6,  -0.5,  -0.5 ],
                                        [ -1.4,  -1.3,  -0.8,  -0.7,  -0.6,  -0.5 ]])

        foa = Akima1DInterpolator((self.A_Conex), self.Offset_A_Conex, axis=1)
        self.Offset_A = foa((self.A_array))
        fob = Akima1DInterpolator((self.A_Conex), self.Offset_B_Conex, axis=1)
        self.Offset_B = fob((self.A_array))
        fta = Akima1DInterpolator((self.A_Conex), self.Tau_A_Conex, axis=1)
        self.Tau_A = fta((self.A_array))
        ftb = Akima1DInterpolator((self.A_Conex), self.Tau_B_Conex, axis=1)
        self.Tau_B = ftb((self.A_array))
        #print 'self.Offset_A.shape', self.Offset_A.shape
        #print 'self.Offset_A.shape', self.Offset_A.shape
        #print 'self.Offset_B.shape', self.Offset_B.shape
        #print 'self.Tau_A.shape',    self.Tau_A.shape
        #print 'self.Tau_B.shape',    self.Tau_B.shape
        # Create Energy Dependet Parameters
        
        self.Offset_A = np.einsum('ij,k->ijk', self.Offset_A,  np.ones(len(self.log10_Energy_array)))
        self.Offset_B = np.einsum('ij,k->ijk', self.Offset_B,  (self.log10_Energy_array-19.))
        self.Offset = self.Offset_A + self.Offset_B
        
        self.Tau_A = np.einsum('ij,k->ijk', self.Tau_A,  np.ones(len(self.log10_Energy_array)))
        self.Tau_B = np.einsum('ij,k->ijk', self.Tau_B,  (self.log10_Energy_array-19.))
        self.Tau = self.Tau_A + self.Tau_B
        # Extend Grid to all nuclear masses from 1 to 56 by interpolating the Conex results
        
        #print 'self.Offset.shape', self.Offset.shape
        #print 'self.Tau.shape',    self.Tau.shape
        #print len(self.A_array), self.A_array
        
        self.Mean_Matrix = self.Offset + 5.*self.Tau
        self.MeanSquare_Matrix  = 30.*self.Tau**2 + self.Offset**2 + 10.*self.Tau*self.Offset
        
        #print 'self.Mean_Matrix.shape', self.Mean_Matrix.shape
        #print 'self.MeanSquare_Matrix.shape', self.MeanSquare_Matrix.shape
        # model grid is 2 x 56 x len(Energies)
        
    def getMeanRMS(self, f_A_array):
        Mean       = np.einsum('ijk,jk->ik', self.Mean_Matrix,       f_A_array)
        MeanSquare = np.einsum('ijk,jk->ik', self.MeanSquare_Matrix, f_A_array)
        #print Mean.shape

        #Mean       = np.einsum('ijk,j', self.Mean_Matrix,       f_A_array)
        #MeanSquare = np.einsum('ijk,j', self.MeanSquare_Matrix, f_A_array)
        RMS = np.sqrt(MeanSquare - Mean**2)
        return Mean, RMS


    
