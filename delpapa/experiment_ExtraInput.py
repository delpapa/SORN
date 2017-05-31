from __future__ import division
from pylab import *
import utils
utils.backup(__file__)

from delpapa.plot import plot_results as plot_results_single

from common.sources import CountingSource, NoSource
from common.experiments import AbstractExperiment
from common.sorn_stats import *

class Experiment_test(AbstractExperiment):
    def start(self):
        super(Experiment_test,self).start()
        c = self.params.c

        self.inputsource = NoSource(N_i=c.N_u_e)

        stats_single = [
                         ActivityStat(),
                         SpikesStat(),
                         ConnectionFractionStat(),
                        ]
        return (self.inputsource, stats_single)

    def reset(self,sorn):
        super(Experiment_test,self).reset(sorn)
        c = self.params.c
        stats = sorn.stats # init sets sorn.stats to None
        sorn.__init__(c,self.inputsource)
        sorn.stats = stats

    def run(self,sorn):
        super(Experiment_test,self).run(sorn)
        c = self.params.c

        print '\n\nTransient...'
        sorn.simulation(c.steps_transient)

        print '\n\nExternal Input off:'
        sorn.simulation(c.steps_noExternalInput)


        print '\n\nExternal Input on'

        # external input definition
        word1 = "A"
        word2 = "B"
        word3 = "C"
        word4 = "D"
        word5 = "E"
        word6 = "F"
        word7 = "G"
        word8 = "H"
        word9 = "I"
        word10 = "J"
        words = [word1, word2, word3, word4, word5, word6, word7, \
                word8, word9, word10]
        # extrenal input random transitions
        words_len = len(words)
        words_trans = np.ones([words_len,words_len])*(1./words_len)

        newsource = CountingSource(words, words_trans, \
                                            c.N_u_e,c.N_u_i,avoid=False)
        sorn.source = newsource
        sorn.W_eu = newsource.generate_connection_e(c.N_e)
        sorn.simulation(c.steps_ExternalInput)

        return {'source_plastic':self.inputsource,'source_test':newsource}

    def plot_single(self,path,filename):
        plot_results_single(path,filename)
