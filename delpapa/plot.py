from __future__ import division
from pylab import *
from scipy.optimize import curve_fit
from scipy import stats
from scipy import signal
import tables

import sys
sys.path.insert(0,"../")
import utils
utils.backup(__file__)

from scipy.io import savemat
import datetime
from utils.pca import pca
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats.stats import pearsonr

import cPickle as pickle
import gzip
from common.sources import TrialSource
import os
import platform

plot_MDS = True
normalize_PCA = False
matplotlib.rcParams.update({'font.size': 10}) 

def plot_results(result_path,result):   
    h5 = tables.openFile(os.path.join(result_path,result),'r')
    data = h5.root 
    pickle_dir = data.c.logfilepath[0]
    plots_path = os.path.join('..','plots')
    if not os.path.exists(plots_path):
        os.mkdir(plots_path)
    os.chdir(plots_path)
        
	### Plot the connection fraction
    if data.__contains__('ConnectionFraction'):
        print 'plot connectionfraction'
        figure()
        plot(data.ConnectionFraction[0][:data.c.N_steps[0]])
        xlabel('Time Step'); ylabel('Fraction of E-E connections')
        tight_layout()
        utils.saveplot('ConnectionFraction.pdf')


    if data.__contains__('activity') and False:
        print 'plot activity'
        figure()
        activity = data.activity[0, :]
        plot(activity, 'k')
        xlabel('Step'); ylabel('activity')
        tight_layout()
        utils.saveplot('Activity.pdf')

	### Plot the E spikes
    ### plot the raster of the last last_n_spikes steps	
    ### plot the activity of the last last_n_spikes steps			 
    if data.__contains__('Spikes') and False:
      
        print 'plot spikes'
        # raster plot (last_n_spikes)
        last_spikes = data.c.stats.only_last_spikes[0]
        spikes = data.Spikes[0, :, -last_spikes:]
        figure()
        if data.__contains__('activity'):
			subplot(211)
        steps = -1 # data.c.steps_plastic[0]
        for (i,sp) in enumerate(spikes):
            s_train = where(sp == 1)[0]
            if s_train != []:
				vlines(s_train, i + 0.5, i + 1.5)
				hold('on')
        ylabel('Excitatory Neuron')
        ylim(0.5, i + 1.5)
        
        if data.__contains__('activity'):
			activity = data.activity[0, -last_spikes:]
			subplot(212)
			plot(activity, 'k')
			xlabel('Step'); ylabel('activity')
        tight_layout()
        utils.saveplot('Raster_end.pdf')
          
    if data.__contains__('SpikesInh'):
		
        print 'plot spikesInh'
        last_spikes = data.c.stats.only_last_spikes[0]       
        spikes = data.SpikesInh[0, :, -last_spikes:]
        figure()
        if data.__contains__('activity'):
			subplot(211)
        steps = -1 # data.c.steps_plastic[0]
        for (i,sp) in enumerate(spikes):
            s_train = where(sp == 1)[0]
            if s_train != []:
				vlines(s_train, i + 0.5, i + 1.5)
				hold('on')
        ylabel('Inhibitory Neuron')
        ylim(0.5, i + 1.5)
        if data.__contains__('activityInh'):
			activity = data.activityInh[0, -last_spikes:]
			subplot(212)
			plot(activity, 'k')
			xlabel('Step'); ylabel('activity')
        tight_layout()
        utils.saveplot('Raster_end_inh.pdf')
    
    ### Plot the degree distribution and the curve fit of the end weight
    if data.__contains__('endweight'):
        print 'plot endweight'
        		
        N_e = data.c.N_e[0]     

        # First the logweight:
        logweight = data.endweight[0][data.endweight[0]>0]
        figure()
        logbins = logspace(-2,0,10)
        (y,_) = histogram(logweight,bins=logbins)
        #fit data to lognormal
        x = logbins[:-1]+(logbins[0]+logbins[1])/2.0
        semilogx(x,y,'.')

        # Do the fitting
        def lognormal(x,mue,var,scale):
            return scale * (exp(- ((log(x)-mue)*(log(x)-mue)) /\
                                          (2*var)) / (x*sqrt(2*pi*var)))

        try:
            popt, pcov = curve_fit(lognormal, x, y)
        except(RuntimeError):
            popt = [0,0,0]
        curve_x = logspace(-2,0,100)
        fitted_y = lognormal(curve_x,*popt)
        plot(curve_x,fitted_y)
        xlabel('Weight'); ylabel('Frequency')
        legend(('Data', 'Logn. fit'))
        tight_layout()
        utils.saveplot('LogWeights_%s.pdf'%\
                                          (data.c.stats.file_suffix[0]))
        
        # Now scale-free property
        tmp = data.endweight[0]>0.
        binary_connections = tmp
        in_degree = sum(binary_connections,1)
        out_degree = sum(binary_connections,0)
        fig = figure()
        fig.add_subplot(131)
        hist(in_degree)
        ylabel('Frequency')
        xlabel('In degree')
        fig.add_subplot(132)
        hist(out_degree)
        xlabel('Out degree')
        fig.add_subplot(133)
        hist(in_degree+out_degree)
        xlabel('In+out degree')
        plt.suptitle('Degree distributions')
        utils.saveplot('Degree_Distributions_%s.pdf'%\
                                          (data.c.stats.file_suffix[0]))
                                          
        # Now the weight matrix
        fig_mat = figure()
        imshow(data.endweight[0], cmap='Greys', interpolation='nearest')
        colorbar()
        tight_layout()
        utils.saveplot('FinalWeightMatrix_%s.pdf'%\
                                          (data.c.stats.file_suffix[0])) 
    
    ### Plot ISIs
    if data.__contains__('ISIs'):
        print 'plot ISIs'
        
        figure()
        ISIs = data.ISIs[0]
        x = np.array(range(0,shape(ISIs)[1]))
        y = ISIs[randint(0,shape(ISIs)[0])]
        

        # Do the fitting
        def exponential(x, a, b):
            return a * np.exp(-b*x)
        if data.c.stats.__contains__('ISI_step'):
            start = int(round(mean(argmax(data.ISIs[0],1))))
            x_fit = x[start::data.c.stats.ISI_step[0]]
            y_fit = y[start::data.c.stats.ISI_step[0]]
        else:
            x_fit = x
            y_fit = y
        popt, pcov = curve_fit(exponential, x_fit, y_fit)
        x = np.array(range(shape(ISIs)[1]))
        fitted_y = exponential(x,*popt)
        if data.c.stats.__contains__('ISI_step'):
            plot(x[start+1:],y[start:], '.')
            plot(x[start:],fitted_y[start:])
            xlim([start-1,max(x)])
        else:
            plot(x,y, '.')
            plot(x,fitted_y)
            xlim([-1,shape(ISIs)[1]])
        #~ x_lim[0] -= 1
        #~ xlim(x_lim)

        #~ title('Interspike Intervals')# (%s)'%(data.c.stats.file_suffix[0]))
        xlabel('ISI')
        ylabel('Frequency')
        if data.c.stats.__contains__('ISI_step'):
            fitlabel = 'Exp. fit (sampling=%d)'%data.c.stats.ISI_step[0]
        else:
            fitlabel = 'Exp. fit'
        legend(('Data', fitlabel))# (scale:%.3f exponent:%.3f)'%(popt[0],-popt[1])))
        tight_layout()
        utils.saveplot('ISIs_%s.pdf'%(data.c.stats.file_suffix[0]))
        
        # Plot the CV hist
        # TODO CHANGE
        #~ if data.c.stats.__contains__('ISI_step'):
            #~ ISIs = ISIs[:,start::data.c.stats.ISI_step[0]]
        figure()
        CVs = []
        intervals = arange(shape(ISIs)[1])
        for i in range(shape(ISIs)[0]):
            tmp = repeat(intervals.tolist(),ISIs[i].tolist())
            CVs.append(std(tmp)/mean(tmp))
        CVs = array(CVs)
        CVs = CVs[CVs>0] #ignore nans
        n, bins, patches = hist(CVs)
        xlim([floor((min(bins[n>0])-0.1)*10)*0.1, 
                                    floor((max(bins[n>0])+0.2)*10)*0.1])
        xlabel('ISI CV')
        ylabel('Frequency')
        utils.saveplot('ISIs_CV_%s.pdf'%(data.c.stats.file_suffix[0]))
    
    
    h5.close()
    tight_layout()
    # show()

def plot_results_perturbation(result_path,result):
    h5 = tables.openFile(os.path.join(result_path,result),'r')
    data = h5.root 
    pickle_dir = data.c.logfilepath[0]
    plots_path = os.path.join('..','plots')
    if not os.path.exists(plots_path):
        os.mkdir(plots_path)
    os.chdir(plots_path)
    	
	### Plot the connection fraction
    if data.__contains__('ConnectionFraction'):
        print 'plot connectionfraction'
        figure()
        non_pert_data = data.ConnectionFraction[0][:data.c.steps_plastic[0]+data.c.steps_perturbation[0]]
        pert_data = data.ConnectionFraction[0][data.c.steps_plastic[0]+data.c.steps_perturbation[0]:]
        
        plot(non_pert_data, label='Without perturbation')
        non_pert_data_x = np.arange(size(non_pert_data)-size(pert_data), size(non_pert_data))
        plot(non_pert_data_x, pert_data, 'r', label='With perturbation')
        xlabel('Time Step'); ylabel('Fraction of E-E connections')
        legend(loc='best')
        tight_layout()
        utils.saveplot('ConnectionFraction.pdf')

    ### plot the activity difference after the perturbation
    if data.__contains__('activity') and False:
        print 'plot activity'
        figure()
        non_pert_act = data.activity[0][:data.c.steps_plastic[0]+data.c.steps_perturbation[0]]
        pert_act = data.activity[0][data.c.steps_plastic[0]+data.c.steps_perturbation[0]:]
        
        plot(non_pert_act, label='Without perturbation')
        non_pert_act_x = np.arange(size(non_pert_act)-size(pert_act), size(non_pert_act))
        plot(non_pert_act_x, pert_act, 'r', label='With perturbation')
        xlabel('Time Step'); ylabel('activity')
        legend(loc='best')
        tight_layout()

	### Plot the E spikes
    ### plot the raster of the non-perturbated and perturbated spikes	
    ### plot the difference of spikes			 
    if data.__contains__('Spikes'):
      
        print 'plot spikes'
        # raster plot (last_n_spikes)
        last_spikes = data.c.stats.only_last_spikes[0]
        non_pert_spikes = data.Spikes[0, :, -last_spikes:-last_spikes/2]
        pert_spikes = data.Spikes[0, :, -last_spikes/2:]
        
        only_non_pert = non_pert_spikes - pert_spikes
        only_non_pert[only_non_pert < 0] = 0
        only_pert = pert_spikes - non_pert_spikes
        only_pert[only_pert < 0] = 0 
              
        figure()
        subplot(211)
        steps = -1 # data.c.steps_plastic[0]
        for (i,sp) in enumerate(only_non_pert):
            s_train = where(sp == 1)[0]
            if s_train != []:
                if i == 1: # small change to add the label - TODO looks terrible, change this ASAP 
                    vlines(s_train, i + 0.5, i + 1.5, 'b', label = 'Without perturbation only')
                else:
                    vlines(s_train, i + 0.5, i + 1.5, 'b')
                hold('on')
        for (i,sp) in enumerate(only_pert):
            s_train = where(sp == 1)[0]
            if s_train != []:
                if i == 1: # small change to add the label - TODO looks terrible, change this ASAP
                    vlines(s_train, i + 0.5, i + 1.5, 'r', label = 'With perturbation only')
                else:
                    vlines(s_train, i + 0.5, i + 1.5, 'r')
                hold('on')
        ylabel('Excitatory Neuron')
        legend(loc='best')
        xlim([0,1000])

        ylim(0.5, i + 1.5)
        
        subplot(212)
        diff = abs(non_pert_spikes-pert_spikes)
        plot(sum(diff, 0), 'k')
        ylabel('# different spikes')
        xlabel('Step after perturbation')
        tight_layout()
        utils.saveplot('Raster_end.pdf')
        
    h5.close()
    tight_layout()
    #~ show() ### this line print the figures on the screen (TURN IF OFF for loop runnings)


#~ if __name__=='__main__':        
    #~ plot_results(r'/home/delpapa/Desktop/sorn/py/backup/test_single/2014-08-27 11-05-36/common','result.h5')
    #~ show()
