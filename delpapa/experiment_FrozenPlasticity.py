from __future__ import division
from pylab import *
import utils
utils.backup(__file__)

from delpapa.plot import plot_results_perturbation as plot_results_single

from common.sources import CountingSource, TrialSource, NoSource
from common.experiments import AbstractExperiment
from common.sorn_stats import *

class Experiment_test(AbstractExperiment):
    def start(self):
        super(Experiment_test,self).start()   
        c = self.params.c
                                  
        self.inputsource = NoSource()
                
        stats_all = []
        stats_single = [
                         ActivityStat(),
                         ConnectionFractionStat(),
                        ]
        return (self.inputsource,stats_all+stats_single,stats_all)
        
    def reset(self,sorn):
        super(Experiment_test,self).reset(sorn)
        c = self.params.c
        stats = sorn.stats # init sets sorn.stats to None
        sorn.__init__(c,self.inputsource)
        sorn.stats = stats
            
    def run(self,sorn):
        super(Experiment_test,self).run(sorn)
        c = self.params.c
        
        
        # Initial pphse with plsaticity
        print '\n\nInitial Phase:'
        sorn.simulation(c.steps_plastic)
        newseed = randint(999999) # random seed 
        tmpstats = sorn.stats
        sorn.stats = 0
        filename = utils.logfilename("net_before_pert.pickle")
        sorn.quicksave(filename)
        sorn.stats = tmpstats
        
        # Reset point: next line runs the 'normal' SORN
        print '\n\nNon-frozen steps:'
        seed(newseed)
        sorn.simulation(c.steps_perturbation)


        # Return to the reset point, freezes plasticity and run again
        print '\n\nFrozen steps:'
        sorn = sorn.quickload(filename)
        sorn.stats = tmpstats
        sorn.stats.obj = sorn

        # Freeze plasticity
        # comment a line NOT to freeze a specific plasticity mechanism
        sorn.W_ee.c.eta_stdp = 0     # freezes STDP
        sorn.W_ei.c.eta_istdp = 0    # freezes iSTDP
        sorn.W_ee.c.sp_prob = 0      # freezes SP
        c.eta_ip = 0                 # freezes IP
        c.noise_sig = 0              # freezes noise
            
        seed(newseed)
        sorn.simulation(c.steps_perturbation)
        
        return {'source_plastic':self.inputsource}
     
    def plot_single(self,path,filename):
        plot_results_single(path,filename)

