from __future__ import division
import numpy as np
import utils
utils.backup(__file__)
from common.defaults import *
import random

########################################################################
# Simple SORN run - for avalanches with the standard parameters        #
# This generates data to reproduce Fig. 1, 2, 3, 5, and S2 of          #
# Del Papa et al. (2017)                                               #
#                                                                      #
# Fig 5: c.noise_sig: (Gaussian noise standart deviation) is used to   #
# vary the noise level. c.noise_fire: if not 0, gives the change of    #
# randomly activating a ramdonly chosen neuron.                        #
#                                                                      #
# Fig. S2: c.noise_fire_struc: if 1, always makes the SAME set of      #
# E neurons fire. Change this number at line 121 of common/sorn.py     #
#                                                                      #  
########################################################################


c.N_e = 200 # number of excitatory neurons 
c.N_i = int(np.floor(0.2*c.N_e)) # percent of inhibitory neurons
c.N = c.N_e + c.N_i # total number of neurons
c.N_u_e = 0 # no input 

c.W_ee = utils.Bunch(use_sparse=True, 
                     lamb = 0.1*c.N_e,
                     avoid_self_connections=True,
                     eta_stdp = 0.004,
                     sp_prob =  c.N_e*(c.N_e-1)*(0.1/(200*199)), #0.1 for N_e = 200, 
                     sp_initial = 0.001, 
                     no_prune = False, 
                     upper_bound = 1
                     )

c.W_ei = utils.Bunch(use_sparse=False,
                     lamb=0.2*c.N_e,
                     avoid_self_connections=True,
                     eta_istdp = 0.001,  
                     h_ip=0.1)  

c.W_ie = utils.Bunch(use_sparse=False,
                     lamb=1.0*c.N_i,
                     avoid_self_connections=True)
                     
c.steps_plastic = 5000000
c.N_steps = c.steps_plastic
c.eta_ip = 0.01
c.h_ip = 0.1
c.noise_sig =  np.sqrt(0.05)  # Phengsheng np.sqrt(0.05) 
c.noise_fire = 0
c.noise_fire_struc = 0 # if 1, activates fixed neurons all time steps                                       


c.stats.file_suffix = 'test'
c.display = True
c.stats.only_last_spikes = 1000 # number of last spikes to save
c.stats.save_spikes = True

c.experiment.module = 'delpapa.experiment_PhengshengSORN'
c.experiment.name = 'Experiment_test'

