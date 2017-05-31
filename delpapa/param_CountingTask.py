from __future__ import division
import numpy as np
import utils
utils.backup(__file__)
from common.defaults import *
import random

c.N_e = 200 # number of excitatory neurons
c.N_i = int(np.floor(0.2*c.N_e)) # percent of inhibitory neurons
c.N = c.N_e + c.N_i # total number of neurons
c.N_u_e = int(np.floor(0.05*c.N_e)) # with input
c.N_u_i = 0

c.W_ee = utils.Bunch(use_sparse=True,
                     lamb = 0.1*c.N_e,
                     avoid_self_connections=True,
                     eta_stdp = 0.004,
                     sp_prob =  0,#########c.N_e*(c.N_e-1)*(0.1/(200*199)), #0.1 for N_e = 200,
                     sp_initial = 0, ### 0.001,
                     no_prune = False,
                     upper_bound = 1
                     )

c.W_ei = utils.Bunch(use_sparse=False,
                     lamb=1*c.N_e,###0.2*c.N_e,
                     avoid_self_connections=True,
                     eta_istdp = 0,###0.001,
                     h_ip=0.1)

c.W_ie = utils.Bunch(use_sparse=False,
                     lamb=1.0*c.N_i,
                     avoid_self_connections=True)

c.steps_plastic = 5000
c.steps_readouttrain = 5000
c.steps_readouttest = 5000
c.steps_avalanches = 0
c.N_steps =  c.steps_plastic + c.steps_readouttrain + \
                               c.steps_readouttest + c.steps_avalanches
c.C_steps = c.steps_readouttrain + c.steps_readouttest

c.eta_ip = 0.01
c.h_ip = 0.1
c.noise_sig = np.sqrt(0.05)
c.noise_fire = 0
c.noise_fire_struc = 0


c.stats.file_suffix = 'test'
c.display = True
c.stats.only_last_spikes = 1000 # number of last spikes to save
c.stats.save_spikes = True
c.input_gain = 1 # needs to be 1!

c.experiment.module = 'delpapa.experiment_CountingTask'
c.experiment.name = 'Experiment_test'
