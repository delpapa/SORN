from __future__ import division
from pylab import *
import utils
utils.backup(__file__)

from delpapa.plot import plot_results as plot_results_single

from common.sources import NoSource
from common.experiments import AbstractExperiment
from common.sorn_stats import *

class Experiment_test(AbstractExperiment):
    def start(self):
        super(Experiment_test,self).start()
        c = self.params.c

        self.inputsource = NoSource()

        stats_single = [
                         ActivityStat(),
                         SpikesStat(),
                         ConnectionFractionStat(),
                        ]

        return (self.inputsource,stats_single)

    def reset(self,sorn):
        super(Experiment_test,self).reset(sorn)
        c = self.params.c
        stats = sorn.stats
        sorn.__init__(c,self.inputsource)
        sorn.stats = stats

    # Run an experiment with plasticity on
    def run(self,sorn):
        super(Experiment_test,self).run(sorn)
        c = self.params.c
        sorn.simulation(c.steps_plastic)

        return {'source_plastic':self.inputsource}

    def plot_single(self,path,filename):
        plot_results_single(path,filename)
