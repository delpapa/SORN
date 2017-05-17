####
# Script for the S2 paper figure (from PhaseTransitionNoise)
# A,B: Fixed spikes 
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
from matplotlib import gridspec # for the different subplot sizes

transient_steps = 2e6
stable_steps = 3e6
number_of_files = 100   
THETA_P = 15           # here, half doesnt work (cause the activity is                                       

### figure parameters
width  =  7
height = 3
fig = figure(1, figsize=(width, height))
fig_S1a = subplot(121)
fig_S1b = subplot(122)

letter_size = 10
letter_size_panel = 12
line_width = 1.5
line_width_fit = 2.0
subplot_letter = (-0.18, 1.05)

########################################################################
# Fig. S1A,B: Fixed Spikes (power-laws and activity)

experiment_name = 'FixedSpikes'
print experiment_name

for p in ['000', '005', '010']:

    print p
    data_all = np.zeros((number_of_files, 3e6)) 
    act_density = zeros((number_of_files, 200)) 
        
    for data_file in xrange(number_of_files):
        exper = 'result.h5'
        exper_path =  '../Avalanche_Experiments/PhaseTransition_Noise/'\
         + experiment_name + '/p' + p + '/' + str(data_file+1)\
                                                            + '/common/'
        h5 = tables.openFile(os.path.join(exper_path,exper),'r')
        data = h5.root
        data_all[data_file] = np.around(data.activity[0] \
                                       [transient_steps:]*data.c.N_e)    
            
        for i in xrange(int(stable_steps)):
            act_density[data_file, data_all[data_file, i]] += 1
        act_density[data_file, :] /= act_density[data_file, :].sum()
            
        h5.close()
            
    # activity distribution and std    
    act_density_mean = act_density.mean(0)
    act_density_std = act_density.std(0) 
            
    # calculates avalanches
    
    # work around to exclude always active units
    if p == '005':
        Act_thresh = int(data_all.mean()/2. + 5)
    elif p == '000':
        Act_thresh = int(data_all.mean()/2.)
    elif p == '010':
        Act_thresh = int(data_all.mean()/2. + 10)    
        
    T_data, S_data = analysis.avalanches(data_all, \
                                       'N', '200', Threshold=Act_thresh)
    subplot(121)
    if p == '005':
        pl.plot_pdf(S_data, linewidth=line_width)
    elif p == '000':
        pl.plot_pdf(S_data, linewidth=line_width)
    elif p == '010':
        pl.plot_pdf(S_data, linewidth=line_width)
        
    xscale('log'); yscale('log')
    xlabel(r'$S$', fontsize=letter_size); 
    ylabel(r'$f(S)$', fontsize=letter_size);
    xlim([1, 3000])
    ylim([0.00001, 0.1])
    xticks([1, 10, 100, 1000], \
        ['$10^0$', '$10^{1}$', '$10^{2}$', '$10^{3}$'])
    yticks([0.1, 0.001, 0.00001],\
            ['$10^{-1}$', '$10^{-3}$', '$10^{-5}$'])
    fig_S1a.spines['right'].set_visible(False)
    fig_S1a.spines['top'].set_visible(False)
    fig_S1a.xaxis.set_ticks_position('bottom')
    fig_S1a.yaxis.set_ticks_position('left')
    tick_params(labelsize=letter_size)  
    
     
        
    subplot(122) 
    if p == '000':
        label_p = '0%'
        plot(act_density_mean, linewidth=line_width, label=label_p)
    elif p == '005':
        label_p = '5%'
        plot(act_density_mean, linewidth=line_width, label=label_p)
    elif p == '010':
        label_p = '10%'
        plot(act_density_mean, linewidth=line_width, label=label_p) 
        
        
    xlabel(r'$a(x)$ [# neurons]', fontsize=letter_size) 
    ylabel(r'$p(a(x))$', fontsize=letter_size)
    xlim([0, 70])
    xticks([0, 20, 40, 60], \
       ['$10$', '$20$', '$40$', '$60$'])
    ylim([0, 0.06])
    yticks([0, 0.05],
       ['0', '0.05'])
    fig_S1b.spines['right'].set_visible(False)
    fig_S1b.spines['top'].set_visible(False)
    fig_S1b.xaxis.set_ticks_position('bottom')
    fig_S1b.yaxis.set_ticks_position('left')
    tick_params(labelsize=letter_size)  
    legend(loc=(0.7, 0.65), prop={'size':letter_size}, \
                                 title=r'Active neurons', frameon=False)
    fig_S1b.get_legend().get_title().set_fontsize(letter_size)
########################################################################
fig_S1a.annotate('A', xy=subplot_letter, xycoords='axes fraction', \
                fontsize=letter_size_panel ,  fontweight='bold', \
                horizontalalignment='right', verticalalignment='bottom') 
fig_S1b.annotate('B', xy=subplot_letter, xycoords='axes fraction', \
                fontsize=letter_size_panel ,  fontweight='bold', \
                horizontalalignment='right', verticalalignment='bottom')
fig.subplots_adjust(wspace=.2)
gcf().subplots_adjust(bottom=0.17)
                 
print 'Saving figures...',		
result_path = '../Avalanche_Results/'
result_name_png = 'FigS3.pdf'
savefig(os.path.join(result_path, result_name_png), format = 'pdf')
print 'done'
