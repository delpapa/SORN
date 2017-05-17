####
# Script for the 5th paper figure (from PhaseTransitionNoise)
# A,B: Gaussian noise 'phase transition'
# C,D: random spike 'phase transition'
# E,F,G: raster plot for every phase 
####

from pylab import * 
import matplotlib.gridspec as gridspec

import tables
import os
from tempfile import TemporaryFile
import scipy, scipy.stats # for the binomial distribution

import data_analysis as analysis

# work around to run powerlaw package [Alstott et al. 2014]
import sys
sys.path.append('/home/delpapa/lib/python2.7/site-packages')
sys.path.append('/home/delpapa/lib/python2.7/site-packages/mpmath')
import powerlaw as pl


transient_steps = 2e6
stable_steps = 3e6
number_of_files = 30
THETA = 'half'
raster_steps = 800

### figure parameters
width  =  6
height = 7
fig = figure(1, figsize=(width, height))
    
fig_5a = subplot(2, 2, 1)
fig_5b = subplot(2, 2, 2)
fig_5c = subplot(2, 2, 3)
fig_5d = subplot(2, 2, 4)

letter_size = 11
letter_size_panel = 13
line_width = 1.5
line_width_fit = 2.0
subplot_letter = (-0.18, 1.05)
subplot_letter_long = (-0.08, 0.85)

########################################################################
# Fig. 5A,B: Gaussian noise (power-laws and activity)

experiment_name = 'Gaussian'
print experiment_name

for sigma in ['0.005', '0.05', '5']:

    print sigma
    data_all = np.zeros((number_of_files, stable_steps)) 
    act_density = zeros((number_of_files, 200)) 
        
    for data_file in xrange(number_of_files):
        exper = 'result.h5'
        exper_path =  '../../../Avalanche_Experiments/PhaseTransition_Noise/'\
         + experiment_name + '/sigma' + sigma + '/' + str(data_file+1)\
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
    T_data, S_data = analysis.avalanches(data_all, \
                                           'N', '200', Threshold=THETA)
    
    if sigma in ['0.005', '0.05']:
                                           
        a_dur1, a_area1 = analysis.avalanches(data_all, 'N', '200',\
                                                      Theta_percent = 5)
        a_dur2, a_area2 = analysis.avalanches(data_all, 'N', '200',\
                                                     Theta_percent = 20)                                         
    
    else:
        a_dur1, a_area1 = analysis.avalanches(data_all, 'N', '200',\
                                                      Threshold=8)
        a_dur2, a_area2 = analysis.avalanches(data_all, 'N', '200',\
                                                     Threshold=11)  
        
                                           
    subplot(221)
    if sigma == '0.05':
        
        a_area1_pdf = pl.pdf(a_area1, 10)
        a_area2_pdf = pl.pdf(a_area2, 10)            
        BinCenters1 = (a_area1_pdf[0][:-1]+a_area1_pdf[0][1:])/2.
        BinCenters2 = (a_area2_pdf[0][:-1]+a_area2_pdf[0][1:])/2.
        x_max = a_area1_pdf[0].max()
        interp1 = np.interp(np.arange(x_max), BinCenters1, a_area1_pdf[1])
        interp2 = np.interp(np.arange(x_max), BinCenters2, a_area2_pdf[1])
        fill_between(np.linspace(0,x_max,len(interp2)), interp2, interp1, facecolor='k', alpha=0.2)
        pl.plot_pdf(S_data, color='k', linewidth=line_width_fit)
        
    elif sigma == '0.005':
        a_area1_pdf = pl.pdf(a_area1, 10)
        a_area2_pdf = pl.pdf(a_area2, 10)            
        BinCenters1 = (a_area1_pdf[0][:-1]+a_area1_pdf[0][1:])/2.
        BinCenters2 = (a_area2_pdf[0][:-1]+a_area2_pdf[0][1:])/2.
        x_max = a_area1_pdf[0].max()
        interp1 = np.interp(np.arange(x_max), BinCenters1, a_area1_pdf[1])
        interp2 = np.interp(np.arange(x_max), BinCenters2, a_area2_pdf[1])
        fill_between(np.linspace(0,x_max,len(interp2)), interp2, interp1, facecolor='darkcyan', alpha=0.2)
        pl.plot_pdf(S_data, color='darkcyan', linewidth=line_width)
        
    elif sigma == '5':
        a_area1_pdf = pl.pdf(a_area1, 10)
        a_area2_pdf = pl.pdf(a_area2, 10)            
        BinCenters1 = (a_area1_pdf[0][:-1]+a_area1_pdf[0][1:])/2.
        BinCenters2 = (a_area2_pdf[0][:-1]+a_area2_pdf[0][1:])/2.
        x_max = a_area1_pdf[0].max()
        interp1 = np.interp(np.arange(x_max), BinCenters1, a_area1_pdf[1])
        interp2 = np.interp(np.arange(x_max), BinCenters2, a_area2_pdf[1])
        fill_between(np.linspace(0,x_max,len(interp2)), interp2, interp1, facecolor='r', alpha=0.2)
        pl.plot_pdf(S_data, color='r', linewidth=line_width)
        
                   
    xscale('log'); yscale('log')
    xlabel(r'$S$', fontsize=letter_size)
    ylabel(r'$f(S)$', fontsize=letter_size)
    # ticks name
    xlim([1, 3000])
    ylim([0.00001, 0.1])
    xticks([1, 10, 100, 1000], \
     ['$10^0$', '$10^{1}$', '$10^{2}$', '$10^{3}$'])
    yticks([0.1, 0.001, 0.00001],\
            ['$10^{-1}$', '$10^{-3}$', '$10^{-5}$'])
    fig_5a.spines['right'].set_visible(False)
    fig_5a.spines['top'].set_visible(False)
    fig_5a.xaxis.set_ticks_position('bottom')
    fig_5a.yaxis.set_ticks_position('left')
    tick_params(labelsize=letter_size)  
            
    subplot(222) 
    if sigma == '0.005':
        plot(act_density_mean, 'c', linewidth=line_width, label='low')
    if sigma == '0.05':
        plot(act_density_mean, 'k', linewidth=line_width_fit, label='intermediate')
    if sigma == '5':
        plot(act_density_mean, 'r', linewidth=line_width, label='high') 
    
    ### Binomial distribution plot
    ### p = 0.1 (\mu_IP); n = 200 (N_E)
    bino_x = np.arange(0, 40)
    bino_y = scipy.stats.binom.pmf(bino_x, 200,0.1)
    subplot(222)
    plot(bino_x, bino_y, '--', color='gray', linewidth=line_width_fit)
    
    xlabel(r'$a(x)$ [# neurons]', fontsize=letter_size)    
    ylabel(r'$p(a(x))$', fontsize=letter_size)
    xlim([0, 70])
    ylim([0, 0.11])
    xticks([0, 20, 40, 60], \
       ['$0$', '$20$', '$40$', '$60$'])
    yticks([0, 0.1], \
       ['$0$', '$0.1$'])
    fig_5b.spines['right'].set_visible(False)
    fig_5b.spines['top'].set_visible(False)
    fig_5b.xaxis.set_ticks_position('bottom')
    fig_5b.yaxis.set_ticks_position('left')
    tick_params(labelsize=letter_size)   
      
    legend(loc=(0.5, 0.55), prop={'size':letter_size}, \
                                 title=r'Noise level', frameon=False)
    fig_5b.get_legend().get_title().set_fontsize(letter_size)              
########################################################################

########################################################################
# Fig. 5C,D: Random Spikes (power-laws and activity)

experiment_name = 'RandomSpikes'
print experiment_name

for p in ['000', '005', '010']:

    print p
    data_all = np.zeros((number_of_files, stable_steps)) 
    act_density = zeros((number_of_files, 200)) 
        
    for data_file in xrange(number_of_files):
        exper = 'result.h5'
        exper_path =  '../../../Avalanche_Experiments/PhaseTransition_Noise/'\
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
            
    # calculates avalanches    
    T_data, S_data = analysis.avalanches(data_all, \
                                           'N', '200', Threshold=THETA)
    
    if p in ['000', '005']:
                                           
        a_dur1, a_area1 = analysis.avalanches(data_all, 'N', '200',\
                                                      Theta_percent = 5)
        a_dur2, a_area2 = analysis.avalanches(data_all, 'N', '200',\
                                                     Theta_percent = 20)                                         
    
    else:
        a_dur1, a_area1 = analysis.avalanches(data_all, 'N', '200',\
                                                      Threshold=9)
        a_dur2, a_area2 = analysis.avalanches(data_all, 'N', '200',\
                                                     Threshold=11)   
    
    subplot(223)
    if p == '005':
        
        a_area1_pdf = pl.pdf(a_area1, 10)
        a_area2_pdf = pl.pdf(a_area2, 10)            
        BinCenters1 = (a_area1_pdf[0][:-1]+a_area1_pdf[0][1:])/2.
        BinCenters2 = (a_area2_pdf[0][:-1]+a_area2_pdf[0][1:])/2.
        x_max = a_area1_pdf[0].max()
        interp1 = np.interp(np.arange(x_max), BinCenters1, a_area1_pdf[1])
        interp2 = np.interp(np.arange(x_max), BinCenters2, a_area2_pdf[1])
        fill_between(np.linspace(0,x_max,len(interp2)), interp2, interp1, facecolor='k', alpha=0.2)
        pl.plot_pdf(S_data, color='k', linewidth=line_width_fit)
    elif p == '000':
        a_area1_pdf = pl.pdf(a_area1, 10)
        a_area2_pdf = pl.pdf(a_area2, 10)            
        BinCenters1 = (a_area1_pdf[0][:-1]+a_area1_pdf[0][1:])/2.
        BinCenters2 = (a_area2_pdf[0][:-1]+a_area2_pdf[0][1:])/2.
        x_max = a_area1_pdf[0].max()
        interp1 = np.interp(np.arange(x_max), BinCenters1, a_area1_pdf[1])
        interp2 = np.interp(np.arange(x_max), BinCenters2, a_area2_pdf[1])
        fill_between(np.linspace(0,x_max,len(interp2))[100:-100], interp2[100:-100], interp1[100:-100], facecolor='darkcyan', alpha=0.2)
        pl.plot_pdf(S_data, color='darkcyan', linewidth=line_width)
        
    elif p == '010':    
        a_area1_pdf = pl.pdf(a_area1, 10)
        a_area2_pdf = pl.pdf(a_area2, 10)            
        BinCenters1 = (a_area1_pdf[0][:-1]+a_area1_pdf[0][1:])/2.
        BinCenters2 = (a_area2_pdf[0][:-1]+a_area2_pdf[0][1:])/2.
        x_max = a_area1_pdf[0].max()
        interp1 = np.interp(np.arange(x_max), BinCenters1, a_area1_pdf[1])
        interp2 = np.interp(np.arange(x_max), BinCenters2, a_area2_pdf[1])
        fill_between(np.linspace(0,x_max,len(interp2)), interp2, interp1, facecolor='r', alpha=0.2)
 
        pl.plot_pdf(S_data, color='r', linewidth=line_width)
    
    xscale('log'); yscale('log')
    xlabel(r'$S$', fontsize=letter_size)
    ylabel(r'$f(S)$', fontsize=letter_size)
    # ticks name
    xlim([1, 3000])
    ylim([0.00001, 0.1])
    xticks([1, 10, 100, 1000], \
     ['$10^0$', '$10^{1}$', '$10^{2}$', '$10^{3}$'])
    yticks([0.1, 0.001, 0.00001],\
            ['$10^{-1}$', '$10^{-3}$', '$10^{-5}$'])
    fig_5c.spines['right'].set_visible(False)
    fig_5c.spines['top'].set_visible(False)
    fig_5c.xaxis.set_ticks_position('bottom')
    fig_5c.yaxis.set_ticks_position('left')
    tick_params(labelsize=letter_size)    
        
    subplot(224) 
    if p == '000':
        label_p = '$0\%$'
        plot(act_density_mean, 'c', linewidth=line_width, label=label_p) 
    elif p == '005':
        label_p = '$5\%$'
        plot(act_density_mean, 'k', linewidth=line_width_fit, label=label_p) 
    if p == '010':
        label_p = '$10\%$'
        plot(act_density_mean, 'r', linewidth=line_width, label=label_p) 
    
    ### Binomial distribution plot
    ### p = 0.1 (\mu_IP); n = 200 (N_E)
    bino_x = np.arange(0, 40)
    bino_y = scipy.stats.binom.pmf(bino_x, 200,0.1)
    subplot(224)
    plot(bino_x, bino_y, '--', color='gray', linewidth=line_width_fit)

    xlabel(r'$a(x)$ [# neurons]', fontsize=letter_size) 
    ylabel(r'$p(a(x))$', fontsize=letter_size)
    xlim([0, 70])
    ylim([0, 0.11])
    xticks([0, 20, 40, 60], \
       ['$0$', '$20$', '$40$', '$60$'])
    yticks([0, 0.1], \
       ['$0$', '$0.1$'])
    fig_5d.spines['right'].set_visible(False)
    fig_5d.spines['top'].set_visible(False)
    fig_5d.xaxis.set_ticks_position('bottom')
    fig_5d.yaxis.set_ticks_position('left')
    tick_params(labelsize=letter_size)       

    # legend stuff
    legend(loc=(0.5, 0.55), prop={'size':letter_size}, \
                                 title=r'$p_{\rm s}$', frameon=False)
    fig_5d.get_legend().get_title().set_fontsize(letter_size)
########################################################################

#~ ######################################################################
#~ # Fig. 5E, F, G: Raster plots
#~ 
#~ subplot(917)
#~ for (i,sp) in enumerate(raster_low):
    #~ s_train = where(sp == 1)[0]
    #~ if s_train != []:
        #~ vlines(s_train, i + 0.5, i + 1.5)
        #~ hold('on')
#~ ylabel('Low', fontsize=letter_size)
#~ ylim([0, 200])
#~ tick_params(axis='both', which='major', labelsize=letter_size) 
#~ xticks([])
#~ yticks([])
#~ 
#~ subplot(918)
#~ for (i,sp) in enumerate(raster_inter):
    #~ s_train = where(sp == 1)[0]
    #~ if s_train != []:
        #~ vlines(s_train, i + 0.5, i + 1.5)
        #~ hold('on')
#~ ylabel('Interm.', fontsize=letter_size)
#~ ylim([0, 200])
#~ tick_params(axis='both', which='major', labelsize=letter_size)
#~ xticks([]) 
#~ yticks([])  
#~ 
#~ subplot(919)
#~ for (i,sp) in enumerate(raster_high):
    #~ s_train = where(sp == 1)[0]
    #~ if s_train != []:
        #~ vlines(s_train, i + 0.5, i + 1.5)
        #~ hold('on')
#~ xlabel('Time step', fontsize=letter_size)
#~ ylabel('High', fontsize=letter_size)
#~ ylim([0, 200]) 
#~ tick_params(axis='both', which='major', labelsize=letter_size) 
#~ yticks([]) 
#~ xticks([0, 400, 800], \
       #~ ['$0$', '$400$', '$800$'])


########################################################################
fig_5a.annotate('A', xy=subplot_letter, xycoords='axes fraction', \
                fontsize=letter_size_panel ,  fontweight='bold', \
                horizontalalignment='right', verticalalignment='bottom') 
fig_5b.annotate('B', xy=subplot_letter, xycoords='axes fraction', \
                fontsize=letter_size_panel ,  fontweight='bold', \
                horizontalalignment='right', verticalalignment='bottom') 
fig_5c.annotate('C', xy=subplot_letter, xycoords='axes fraction', \
                fontsize=letter_size_panel ,  fontweight='bold', \
                horizontalalignment='right', verticalalignment='bottom') 
fig_5d.annotate('D', xy=subplot_letter, xycoords='axes fraction', \
                fontsize=letter_size_panel ,  fontweight='bold', \
                horizontalalignment='right', verticalalignment='bottom') 
#~ fig_5e.annotate('E', xy=subplot_letter_long, xycoords='axes fraction', \
                #~ fontsize=letter_size_panel ,  fontweight='bold', \
                #~ horizontalalignment='right', verticalalignment='bottom') 
#~ fig_5f.annotate('F', xy=subplot_letter_long, xycoords='axes fraction', \
                #~ fontsize=letter_size_panel ,  fontweight='bold', \
                #~ horizontalalignment='right', verticalalignment='bottom')
#~ fig_5g.annotate('G', xy=subplot_letter_long, xycoords='axes fraction', \
                #~ fontsize=letter_size_panel ,  fontweight='bold', \
                #~ horizontalalignment='right', verticalalignment='bottom')                 
fig.subplots_adjust(wspace=.25)  
fig.subplots_adjust(hspace=.45)

print 'Saving figures...',		
result_path = '../../../Avalanche_Results/'
result_name_png = 'Fig5_1.pdf'
savefig(os.path.join(result_path, result_name_png), format = 'pdf')
print 'done\n\n'
