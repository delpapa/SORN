####
# Script for the S1 paper figure (from FrozenPlasticity)
# A, B: iSTDPFrozen and AllbutiSTDPFrozen
# C, D: STDP+SPFrozen and iSTDP+IPFrozen
####

from pylab import * 

import tables
import os
from tempfile import TemporaryFile

import data_analysis as analysis

# work around to run powerlaw package [Alstott et al. 2014]
import sys
sys.path.append('/home/delpapa/lib/python2.7/site-packages')
sys.path.append('/home/delpapa/lib/python2.7/site-packages/mpmath')
import powerlaw as pl

transient_steps = 2e6
frozen_steps = 2e6
number_of_files = 37 
THETA = 'half'                                                

width  =  7
height = 7 # fix this to the golden ratio
fig_4 = figure(1, figsize=(width, height))
fig_4a = subplot(221)
fig_4b = subplot(222)
fig_4c = subplot(223)
fig_4d = subplot(224)

letter_size = 10
letter_size_panel = 12
line_width = 1.5
line_width_fit = 2.0
subplot_letter = (-0.2, 1)

for experiment_name in ['iSTDPFrozen', 'AllButiSTDPFrozen']:
    
    print experiment_name
    result_path ='../Avalanche_Results/Frozen_Plasticity/' + \
                                                   experiment_name + '/'
                                                   
    for regime in ['non-frozen', 'frozen']:
        
        print '\n', regime
        data_all = np.zeros((number_of_files, frozen_steps))    
        
        for result_file in range(number_of_files):
            exper = 'result.h5'
            exper_path =  '../Avalanche_Experiments/Frozen_Plasticity/'\
               + experiment_name + '/' + str(result_file+1) + '/common/'
            h5 = tables.openFile(os.path.join(exper_path,exper),'r')
            data = h5.root
        
            if regime == 'non-frozen':
                data_all[result_file] = np.around(data.activity[0] \
              [transient_steps:frozen_steps+transient_steps]*data.c.N_e)
        
            if regime == 'frozen':
                data_all[result_file] = np.around(data.activity[0] \
                             [frozen_steps+transient_steps:]*data.c.N_e)            
        
            h5.close()
            
        a_dur, a_area = analysis.avalanches(data_all, 'N', '200', \
                                                       Threshold=THETA)


        if experiment_name == 'iSTDPFrozen':

            print 'Fig. S1A, B...\n'

            subplot(221)
            if regime == 'non-frozen':
                pl.plot_pdf(a_dur, color='k', linewidth=line_width_fit)
            if regime == 'frozen':
                pl.plot_pdf(a_dur, color='b', linewidth=line_width)
            subplot(222)
            if regime == 'non-frozen':
                pl.plot_pdf(a_area, color='k', \
                                linewidth=line_width_fit, label ='SORN')
            if regime == 'frozen':
                label_exp = 'Frozen (iSTDP)'              
                pl.plot_pdf(a_area, linewidth=line_width, color='b', \
                                                        label=label_exp)    
    

        if experiment_name == 'AllButiSTDPFrozen':

            subplot(221)
            if regime == 'frozen':
                pl.plot_pdf(a_dur, color='c', linewidth=line_width)
            subplot(222)
            if regime == 'frozen':
                label_exp = 'Frozen (All but iSTDP)'              
                pl.plot_pdf(a_area, color='c', linewidth=line_width, \
                                                        label=label_exp)

    
# axis stuff    
subplot(221)
xscale('log'); yscale('log')
xlim([1, 3000])
ylim([0.0000001, 1])
xlabel(r'$T$', fontsize=letter_size)
ylabel(r'$f(T)$', fontsize=letter_size)
fig_4a.xaxis.set_ticks_position('bottom')
fig_4a.yaxis.set_ticks_position('left')
tick_params(labelsize=letter_size)  
# ticks name
xticks([1, 10, 100, 1000], \
       ['$10^{0}$', '$10^{1}$', '$10^{2}$', '$10^{3}$'])
yticks([1, 0.001, 0.000001], \
       ['$10^{0}$', '$10^{-3}$', '$10^{-6}$'])



subplot(222)
xscale('log'); yscale('log')
xlim([1, 20000])
ylim([0.00000001, 1]) 	 
xlabel(r'$S$', fontsize=letter_size) 
ylabel(r'$f(S)$', fontsize=letter_size)
fig_4b.xaxis.set_ticks_position('bottom')
fig_4b.yaxis.set_ticks_position('left')
tick_params(labelsize=letter_size)  
# ticks name
xticks([1, 100, 10000], \
       ['$10^{0}$', '$10^{2}$', '$10^{4}$'])
yticks([1, 0.001, 0.000001], \
       ['$10^{0}$', '$10^{-3}$', '$10^{-6}$'])

fig_4a.spines['right'].set_visible(False)
fig_4a.spines['top'].set_visible(False)
fig_4b.spines['right'].set_visible(False)
fig_4b.spines['top'].set_visible(False)

legend(loc=(0, 0), prop={'size':letter_size}, frameon=False)	
########################################################################
number_of_files = 36 
print 'Fig. S1C, D...\n'
for experiment_name in ['STDPandSPFrozen', 'iSTDPandIPFrozen']:
    
    print experiment_name
    result_path ='../Avalanche_Results/Frozen_Plasticity/' + \
                                                   experiment_name + '/'
                                                   
    for regime in ['non-frozen', 'frozen']:
        
        print '\n', regime
        data_all = np.zeros((number_of_files, frozen_steps))    
        
        for result_file in range(number_of_files):
            exper = 'result.h5'
            exper_path =  '../Avalanche_Experiments/Frozen_Plasticity/'\
               + experiment_name + '/' + str(result_file+1) + '/common/'
            h5 = tables.openFile(os.path.join(exper_path,exper),'r')
            data = h5.root
        
            if regime == 'non-frozen':
                data_all[result_file] = np.around(data.activity[0] \
              [transient_steps:frozen_steps+transient_steps]*data.c.N_e)
        
            if regime == 'frozen':
                data_all[result_file] = np.around(data.activity[0] \
                             [frozen_steps+transient_steps:]*data.c.N_e)            
        
            h5.close()
            
        a_dur, a_area = analysis.avalanches(data_all, 'N', '200', \
                                                       Threshold=THETA)

        if experiment_name == 'STDPandSPFrozen':

            subplot(223)
            if regime == 'non-frozen':
                pl.plot_pdf(a_dur, color='k', linewidth=line_width_fit)
            if regime == 'frozen':
                pl.plot_pdf(a_dur, color='b', linewidth=line_width)
            subplot(224)
            if regime == 'non-frozen':
                pl.plot_pdf(a_area, color='k', \
                                linewidth=line_width_fit, label ='SORN')
            if regime == 'frozen':
                label_exp = 'Frozen (STDP + SP)'              
                pl.plot_pdf(a_area, linewidth=line_width, color='b', \
                                                        label=label_exp)    
    

        if experiment_name == 'iSTDPandIPFrozen':

            subplot(223)
            if regime == 'frozen':
                pl.plot_pdf(a_dur, color='c', linewidth=line_width)
            subplot(224)
            if regime == 'frozen':
                label_exp = 'Frozen (iSTDP + IP)'              
                pl.plot_pdf(a_area, color='c', linewidth=line_width, \
                                                        label=label_exp)

    
# axis stuff    
subplot(223)
xscale('log'); yscale('log')
xlim([1, 3000])
ylim([0.0000001, 1])
xlabel(r'$T$', fontsize=letter_size)
ylabel(r'$f(T)$', fontsize=letter_size)
fig_4c.xaxis.set_ticks_position('bottom')
fig_4c.yaxis.set_ticks_position('left')
tick_params(labelsize=letter_size)  

# ticks name
xticks([1, 10, 100, 1000], \
       ['$10^{0}$', '$10^{1}$', '$10^{2}$', '$10^{3}$'])
yticks([1, 0.001, 0.000001], \
       ['$10^{0}$', '$10^{-3}$', '$10^{-6}$'])



subplot(224)
xscale('log'); yscale('log')
xlim([1, 20000])
ylim([0.00000001, 1]) 	 
xlabel(r'$S$', fontsize=letter_size) 
ylabel(r'$f(S)$', fontsize=letter_size)
fig_4d.xaxis.set_ticks_position('bottom')
fig_4d.yaxis.set_ticks_position('left')
tick_params(labelsize=letter_size)  
# ticks name
xticks([1, 100, 10000], \
       ['$10^{0}$', '$10^{2}$', '$10^{4}$'])
yticks([1, 0.001, 0.000001], \
       ['$10^{0}$', '$10^{-3}$', '$10^{-6}$'])

fig_4c.spines['right'].set_visible(False)
fig_4c.spines['top'].set_visible(False)
fig_4d.spines['right'].set_visible(False)
fig_4d.spines['top'].set_visible(False)

legend(loc=(0, 0), prop={'size':letter_size}, frameon=False)
########################################################################

fig_4a.annotate('A', xy=subplot_letter, xycoords='axes fraction', \
                fontsize=letter_size_panel ,  fontweight='bold', \
                horizontalalignment='right', verticalalignment='bottom') 
fig_4b.annotate('B', xy=subplot_letter, xycoords='axes fraction', \
                fontsize=letter_size_panel ,  fontweight='bold', \
                horizontalalignment='right', verticalalignment='bottom')
fig_4c.annotate('C', xy=subplot_letter, xycoords='axes fraction', \
                fontsize=letter_size_panel ,  fontweight='bold', \
                horizontalalignment='right', verticalalignment='bottom') 
fig_4d.annotate('D', xy=subplot_letter, xycoords='axes fraction', \
                fontsize=letter_size_panel ,  fontweight='bold', \
                horizontalalignment='right', verticalalignment='bottom')  
                
gcf().subplots_adjust(bottom=0.17)
fig_4.subplots_adjust(wspace=.4) 
    
print 'Saving figures...',		
result_path = '../Avalanche_Results/'
result_name_png = 'FigS2.pdf'
savefig(os.path.join(result_path, result_name_png), format = 'pdf')
print 'done\n\n'
