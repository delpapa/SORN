from __future__ import division
from pylab import *
from scipy import linalg
import scipy.sparse as sp
import scipy.stats as st

import cPickle as pickle
import gzip

import utils
utils.backup(__file__)
from stats import StatsCollection
from synapses import create_matrix

#~ This contains the 




# Intrinsic Plasticity
def ip(T,x,c):
    """
    Performs intrinsic plasticity
    
    Parameters:
        T: array
            The current thresholds
        x: array
            The state of the network
        c: Bunch
            The parameter bunch
    """
    # T += c.eta_ip*(x.sum()/float(c.N_e)-c.h_ip)
    T += c.eta_ip*(x-c.h_ip)
    return T

class Sorn(object):
    """
    Self-Organizing Recurrent Neural Network
    """
    def __init__(self,c,source):
        """
        Initializes the variables of SORN
        
        Parameters:
            c: bunch
                The bunch of parameters
            source: Source
                The input source
        """
        self.c = c
        self.source = source

        # Initialize weight matrices
        # W_to_from (W_ie = from excitatory to inhibitory)
        self.W_ie = create_matrix((c.N_i,c.N_e),c.W_ie)
        self.W_ei = create_matrix((c.N_e,c.N_i),c.W_ei)
        self.W_ee = create_matrix((c.N_e,c.N_e),c.W_ee)
        self.W_eu = self.source.generate_connection_e(c.N_e)
        self.W_iu = self.source.generate_connection_i(c.N_i)

        # Initialize the activation of neurons
        self.x = rand(c.N_e)<c.h_ip
        self.y = zeros(c.N_i) # CHANGE_A was rand(c.N_i)<mean(c.h_ip)
        self.u = source.next()

        # Initialize the pre-threshold variables
        self.R_x = zeros(c.N_e)
        self.R_y = zeros(c.N_i)

        # Initialize thresholds
        if c.ordered_thresholds: # Andreeas version
            self.T_i = (arange(c.N_i)+0.5)*((c.T_i_max-c.T_i_min)/(1.*c.N_i))+c.T_i_min## CHANGE_A was c.T_i_min + rand(c.N_i)*(c.T_i_max-c.T_i_min)
            self.T_e = (arange(c.N_e)+0.5)*((c.T_e_max-c.T_e_min)/(1.*c.N_e))+c.T_e_min# CHANGE_A was c.T_e_min + rand(c.N_e)*(c.T_e_max-c.T_e_min)
            shuffle(self.T_e)         # CHANGE_A: only T_e are shuffled
        else:
            self.T_i = c.T_i_min + rand(c.N_i)*(c.T_i_max-c.T_i_min)
            
            #~ maximum thresholds for all neurons
            #~ self.T_e = ones(c.N_e)*c.T_e_max
            #~ minimum thresholds for all neurons
            #~ self.T_e = ones(c.N_e)*c.T_e_min   
                   
            self.T_e = c.T_e_min + rand(c.N_e)*(c.T_e_max-c.T_e_min)

        # Activate plasticity mechanisms
        self.update = True
        self.stats = None

    def step(self,u_new):
        """
        Performs a one-step update of the SORN
        
        Parameters:
            u_new: array
                The input for this step
        """
        c = self.c
        # Compute new state
        self.R_x = self.W_ee*self.x-self.W_ei*self.y-self.T_e
        if not c.noise_sig == 0:
            self.R_x += c.noise_sig*np.random.randn(c.N_e)
        
        if not c.ff_inhibition_broad == 0:
            self.R_x -= c.ff_inhibition_broad  
             
        x_temp = self.R_x+c.input_gain*(self.W_eu*u_new)

        if c.k_winner_take_all:
            expected = int(round(c.N_e * c.h_ip))
            ind = argsort(x_temp)
            # the next line fails when c.h_ip == 1
            x_new = (x_temp > x_temp[ind[-expected-1]])+0
        else:	
            x_new = (x_temp >= 0.0)+0
            
        # New noise - prob. of each neuron being active 
        if not c.noise_fire == 0:
            x_new += ((np.random.random(c.N_e) < c.noise_fire)+0)
            x_new[x_new > 1] = 1
        # 'Fixed' noise - same neurons are always active
        if c.noise_fire_struc:
            x_new[12:12] = 1 # change number of E neurons for Fig. S2

        if self.c.fast_inhibit:
            x_used = x_new
        else:
            x_used = self.x

        self.R_y = self.W_ie*x_used - self.T_i
        if self.c.ff_inhibition:
            self.R_y += self.W_iu*u_new
        if not c.noise_sig == 0:
            self.R_y += c.noise_sig*np.random.randn(c.N_i)            
        y_new = (self.R_y >= 0.0)+0

        # Apply plasticity mechanisms
        # Always apply IP
        ip(self.T_e,x_new,self.c)
        # Apply the rest only when update==true
        if self.update:
            assert self.sane_before_update()
            self.W_ee.stdp(self.x,x_new)
            self.W_eu.stdp(self.u,u_new,to_old=self.x,to_new=x_new)

            self.W_ee.struct_p()
            self.W_ei.istdp(self.y,x_new)
            
            self.synaptic_scaling()

            assert self.sane_after_update()

        self.x = x_new
        self.y = y_new
        self.u = u_new
        
        # Update statistics
        self.stats.add()
    
    def synaptic_scaling(self):
        """
        Performs synaptic scaling for all matrices
        """
        self.W_ee.ss()
        self.W_ei.ss() # this was also found in the EM study
        if self.W_eu.c.has_key('eta_stdp') and self.W_eu.c.eta_stdp>0:
            self.W_eu.ss()

    def sane_before_update(self):
        """
        Basic sanity checks for thresholds and states before plasticity
        """
        eps = 1e-6
        assert all(isfinite(self.T_e))
        assert all(isfinite(self.T_i))

        assert all((self.x==0) | (self.x==1))
        assert all((self.y==0) | (self.y==1))

        # TODO quick and easy fix --> this should be 0.0 but threw
        # assertion sometimes on the cluster
        if self.c.noise_sig == -1.0:
            assert all( self.R_x+self.T_e <= +1.0 + eps)
            assert all( self.R_x+self.T_e >= -1.0 - eps)
            assert all( self.R_y+self.T_i <= +1.0 + eps)
            assert all( self.R_y+self.T_i >= 0.0 - eps)

        return True

    def sane_after_update(self):
        """
        Basic sanity checks for matrices after plasticity
        """
        assert self.W_ee.sane_after_update()
        #assert self.W_ei.sane_after_update()
        assert self.W_ie.sane_after_update()

        return True

    def simulation(self,N,toReturn=[]):
        """
        Simulates SORN for a defined number of steps
        
        Parameters:
            N: int
                Simulation steps
            toReturn: list
                Tracking variables to return. Options are: 'X','Y','R_x'
                'R_y', 'U'
        """
        c = self.c
        source = self.source
        ans = {}

        # Initialize tracking variables
        if 'X' in toReturn:
            ans['X'] = zeros( (N,c.N_e) )
        if 'Y' in toReturn:
            ans['Y'] = zeros( (N,c.N_i) )
        if 'R_x' in toReturn:
            ans['R_x'] = zeros( (N,c.N_e) )
        if 'R_y' in toReturn:
            ans['R_y'] = zeros( (N,c.N_i) )
        if 'U' in toReturn:
            ans['U'] = zeros( (N,source.global_range()) )

        # Simulation loop
        for n in range(N):
            
            # Simulation step
            self.step(source.next())

            # Tracking
            if 'X' in toReturn:
                ans['X'][n,:] = self.x
            if 'Y' in toReturn:
                ans['Y'][n,:] = self.y
            if 'R_x' in toReturn:
                ans['R_x'][n,:] = self.R_x
            if 'R_y' in toReturn:
                ans['R_y'][n,:] = self.R_y
            if 'U' in toReturn:
                ans['U'][n,source.global_index()] = 1

            # Command line progress message
            if c.display and (N>100) and \
                    ((n%((N-1)//100) == 0) or (n == N-1)):
                sys.stdout.write('\rSimulation: %3d%%'%((int)(n/(N-1)\
                                                                *100)))
                sys.stdout.flush()
        return ans

    def quicksave(self, filename=None):
        """
        Saves this object
        
        Parameters:
            filename: string
                Filename to save in. Default: "net.pickle"
        """
        if filename == None:
            filename = utils.logfilename("net.pickle")
        temp = self.stats
        pickle.dump( self, gzip.open( filename, "wb" ), \
                     pickle.HIGHEST_PROTOCOL )
        self.stats = temp

    @classmethod
    def quickload(cls, filename):
        """
        Loads a SORN
        
        Parameters:
            filename: string
                File to load from.
        """
        return pickle.load(gzip.open(filename, 'rb'))


