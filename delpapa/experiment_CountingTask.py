from __future__ import division
from pylab import *
import utils
utils.backup(__file__)

from delpapa.plot import plot_results as plot_results_single

from common.sources import CountingSource, TrialSource, NoSource
from common.experiments import AbstractExperiment
from common.sorn_stats import *

n_middle = 4

class Experiment_test(AbstractExperiment):
    def start(self):
        super(Experiment_test,self).start()
        c = self.params.c

        word1 = 'A'
        word2 = 'D'
        for i in range(n_middle):
            word1 += 'B'
            word2 += 'E'
        word1 += 'C'
        word2 += 'F'
        m_trans = np.ones((2,2))*0.5

        self.inputsource = CountingSource([word1, word2],m_trans,
                           c.N_u_e,c.N_u_i,avoid=False)

        stats_single = [
                         ActivityStat(),
                         CountingLetterStat(),
                         CountingActivityStat(),
                         ConnectionFractionStat(),
                         SpikesStat()
                        ]
        return (self.inputsource,stats_single)

    def reset(self,sorn):
        super(Experiment_test,self).reset(sorn)
        c = self.params.c
        stats = sorn.stats # init sets sorn.stats to None
        sorn.__init__(c,self.inputsource)
        sorn.stats = stats

    def run(self,sorn):
        super(Experiment_test,self).run(sorn)
        c = self.params.c

        #----- Input with plasticity
        print '\nInput plastic period:'
        #~ sorn.W_eu = newsource.generate_connection_e(c.N_e)
        sorn.simulation(c.steps_plastic)

        #----- Input without plasticity - train


        print '\nInput training period:'
        # Turn off plasticity
        sorn.W_ee.c.eta_stdp = 0
        sorn.W_ei.c.eta_istdp = 0
        sorn.W_ee.c.sp_prob = 0
        c.eta_ip = 0
        # turn off noise
        c.noise_sig = 0



        sorn.simulation(c.steps_readouttrain)

        #----- Input without plasticity - test performance
        print '\nInput test period:'
        sorn.simulation(c.steps_readouttest)

        print '\nAvalanche measurement period:'

        # Turn on plasticity
        sorn.W_ee.c.eta_stdp = 0.004
        sorn.W_ei.c.eta_istdp = 0.001
        sorn.W_ee.c.sp_prob = c.N_e*(c.N_e-1)*(0.1/(200*199))
        c.eta_ip = 0.01
        # turn on noise
        c.noise_sig = 0#np.sqrt(0.05)

        sorn.simulation(c.steps_avalanches)

        return {'source_plastic':self.inputsource}

    def plot_single(self,path,filename):
        plot_results_single(path,filename)
