from __future__ import division
import numpy as np
import utils
utils.backup(__file__)
from common.defaults import *

########################################################################
# Perturberd SORN - experiment for the SORN with frozen plasticity     #
# for Fig 4 and S1. This runs the sorn for c.steps_plastic steps,      #
# creates a reset point, run the 'normal' simulation for               #
# c.steps_perturbation steps, and finally resets the parameters back   #
# to the reset point, before running the 'frozen' SORN with frozen     #
# plasticity mechanisms. In order to freeze different plasticity       #
# mechanisms, go to line 61 at experiment_Perturbation.py              #                                                                                                                         #  
########################################################################

c.N_e = 200 # number of excitatory neurons
c.N_i = int(np.floor(0.2*c.N_e)) # percent of inhibitory neurons
c.N = c.N_e + c.N_i # total number of neurons
c.N_u_e = 0 # no input 

c.W_ee = utils.Bunch(use_sparse=True, 
                     lamb=0.1*c.N_e,
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
                     
c.steps_plastic = 2000000 # steps before the perturbation
c.steps_perturbation = 2000000 # steps after the perturbation
c.N_steps = (c.steps_plastic + 2*c.steps_perturbation)
c.N_iterations = 1
c.eta_ip = 0.01
c.h_ip = 0.1
c.noise_sig = np.sqrt(0.05) 
c.noise_fire = 0 # 0.05
c.noise_fire_struc = 0 # 1


c.stats.file_suffix = 'test'
c.display = True
# save the spikes of the perturbation perdiod:
# first half is the non-perturbated network
# second half is the perturbated network
c.stats.only_last_spikes = (2*c.steps_perturbation)
c.stats.save_spikes = True

c.experiment.module = 'delpapa.experiment_FrozenPlasticity'
c.experiment.name = 'Experiment_test'
