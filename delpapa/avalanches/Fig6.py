####
# Script for the 6th paper figure (from ExtraInput)
# A,B: Comparison among SORN, input onset and readaptation
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

section_steps = 2e6
extrainput_steps = 20 # 30 is too much!

number_of_files = 50 #50
number_of_files_trans = 250 #250 # well, this is a lot!
experiment_folder = 'NU(0.02)PlasticityBigGain'
## 'NU(0.02)PlasticityBigGain' - From now on, use this!

possible_regimes = ['normal', 'extrainput_start', 'extrainput_end']

### figure parameters
width  =  7
height = 3
fig_6 = figure(1, figsize=(width, height))

letter_size = 10
letter_size_panel = 12
line_width = 1.5
line_width_fit = 2.0
subplot_letter = (-0.15, 0.9)

########################################################################
# Fig. 1A,B: Gaussian Noise + ABCD + Gaussian Noise

for regime in possible_regimes:
    
    print '\n', regime, '...'
    
    
    
    if regime == 'extrainput_start':
        n_trials = number_of_files_trans
        data_all = zeros((n_trials, extrainput_steps))
        data_all2 = zeros((n_trials, 10*extrainput_steps))
    elif regime == 'extrainput_end':
        n_trials = number_of_files
        data_all = zeros((n_trials, section_steps-extrainput_steps))
    elif regime == 'normal':
        n_trials = number_of_files
        data_all = zeros((n_trials, section_steps))
    
    
    for result_file in xrange(n_trials):

        exper = 'result.h5'
        exper_path =  '../../../Avalanche_Experiments/Extra_Input/' + \
                                                   experiment_folder + \
                                   '/' + str(result_file+1) + '/common/'
        h5 = tables.openFile(os.path.join(exper_path,exper),'r')
        data = h5.root
        
        if regime == 'normal':
            data_all[result_file] = np.around(data.activity[0] \
                             [section_steps:2*section_steps]*data.c.N_e)  
                                              
        if regime == 'extrainput_start':
            # 100 * etra input steps here is to count all avalanches 
            # that start in the window I want
            data_all[result_file] = np.around(data.activity[0] \
                   [2*section_steps:2*section_steps + extrainput_steps]\
                                                            *data.c.N_e)
            data_all2[result_file] = np.around(data.activity[0] \
                [2*section_steps:2*section_steps + 10*extrainput_steps]\
                                                            *data.c.N_e)    
                                                               
        if regime == 'extrainput_end':
            data_all[result_file] = np.around(data.activity[0] \
                       [2*section_steps + extrainput_steps:]*data.c.N_e)                                                                                                             
        h5.close()

    ### Fig3 : Plot activity distribution
    figure(3)
    act_density = zeros((number_of_files, data_all.max()+1))    
    for data_file in xrange(number_of_files):
        
        steps = section_steps
        if regime == 'extrainput_start':
            steps = extrainput_steps
        if regime == 'extrainput_end':
            steps = section_steps-extrainput_steps
        for i in xrange(int(steps)):
            act_density[data_file, data_all[data_file, i]] += 1
        act_density[data_file, :] /= act_density[data_file, :].sum()

    act_density_mean = act_density.mean(0)
    if regime == 'normal':
        plot(act_density_mean, 'b')
    if regime == 'extrainput_start':
        plot(act_density_mean, 'r')
    if regime == 'extrainput_end':
        plot(act_density_mean, 'g')    

    print 'Mean =', data_all.mean()
    print 'Std =', data_all.std()
               
    if regime == 'normal':    
        Thres_normal = int(data_all.mean()/2.) + 1 # rounding purposes
        T_data, S_data = analysis.avalanches(data_all, 'N', '200', \
                                                 Threshold=Thres_normal)
                                                 
        a_dur1, a_area1 = analysis.avalanches(data_all, 'N', '200',\
                                                      Theta_percent = 5)
        a_dur2, a_area2 = analysis.avalanches(data_all, 'N', '200',\
                                                     Theta_percent = 25)                                             
                                                  
    if regime == 'extrainput_start':    
        Thres_start = Thres_normal # 20
        T_data, S_data = analysis.avalanches(data_all2, 'N', '200', \
                    Threshold=Thres_start, Transient = extrainput_steps)
                    
        a_dur1, a_area1 = analysis.avalanches(data_all2, 'N', '200',\
                         Threshold=9, Transient = extrainput_steps)
        a_dur2, a_area2 = analysis.avalanches(data_all2, 'N', '200',\
                          Threshold=12, Transient = extrainput_steps)
    
                                                   
    if regime == 'extrainput_end':  
        Thres_end = Thres_normal
        T_data, S_data = analysis.avalanches(data_all, 'N', '200', \
                                                    Threshold=Thres_end)
        a_dur1, a_area1 = analysis.avalanches(data_all, 'N', '200',\
                                                      Theta_percent = 5)
        a_dur2, a_area2 = analysis.avalanches(data_all, 'N', '200',\
                                                     Theta_percent = 25)                                                
                                                        
                                                                                                
    figure(1)
    fig_6a = subplot(121)
    if regime == 'normal':        
        pl.plot_pdf(T_data, color='k', linewidth = line_width_fit)
        a_dur1_pdf = pl.pdf(a_dur1, 10)
        a_dur2_pdf = pl.pdf(a_dur2, 10) 
        BinCenters1 = (a_dur1_pdf[0][:-1]+a_dur1_pdf[0][1:])/2.
        BinCenters2 = (a_dur2_pdf[0][:-1]+a_dur2_pdf[0][1:])/2.
        x_max = a_dur1_pdf[0].max()
        interp1 = np.interp(np.arange(x_max), BinCenters1, a_dur1_pdf[1])
        interp2 = np.interp(np.arange(x_max), BinCenters2, a_dur2_pdf[1])
        fill_between(np.linspace(0,x_max,len(interp2)), interp2, interp1, facecolor='k', alpha=0.2)

        
        
    if regime == 'extrainput_start':
        a_dur1_pdf = pl.pdf(a_dur1, 10)
        a_dur2_pdf = pl.pdf(a_dur2, 10) 
        BinCenters1 = (a_dur1_pdf[0][:-1]+a_dur1_pdf[0][1:])/2.
        BinCenters2 = (a_dur2_pdf[0][:-1]+a_dur2_pdf[0][1:])/2.
        x_max = a_dur1_pdf[0].max()
        interp1 = np.interp(np.arange(x_max), BinCenters1, a_dur1_pdf[1])
        interp2 = np.interp(np.arange(x_max), BinCenters2, a_dur2_pdf[1])
        fill_between(np.linspace(0,x_max,len(interp2)), interp2, interp1, facecolor='r', alpha=0.2)

        pl.plot_pdf(T_data, color='r', linewidth = line_width)         
    if regime == 'extrainput_end':
        
        a_dur1_pdf = pl.pdf(a_dur1, 10)
        a_dur2_pdf = pl.pdf(a_dur2, 10)            
        BinCenters1 = (a_dur1_pdf[0][:-1]+a_dur1_pdf[0][1:])/2.
        BinCenters2 = (a_dur2_pdf[0][:-1]+a_dur2_pdf[0][1:])/2.
        x_max = a_dur1_pdf[0].max()
        interp1 = np.interp(np.arange(x_max), BinCenters1, a_dur1_pdf[1])
        interp2 = np.interp(np.arange(x_max), BinCenters2, a_dur2_pdf[1])
        fill_between(np.linspace(0,x_max,len(interp2)), interp2, interp1, facecolor='cyan', alpha=0.2)
  
        pl.plot_pdf(T_data, color='cyan', linewidth = line_width_fit)

    fig_6b = subplot(122)
    if regime == 'normal':
        
        a_area1_pdf = pl.pdf(a_area1, 10)
        a_area2_pdf = pl.pdf(a_area2, 10) 
        BinCenters1 = (a_area1_pdf[0][:-1]+a_area1_pdf[0][1:])/2.
        BinCenters2 = (a_area2_pdf[0][:-1]+a_area2_pdf[0][1:])/2.
        x_max = a_area1_pdf[0].max()
        interp1 = np.interp(np.arange(x_max), BinCenters1, a_area1_pdf[1])
        interp2 = np.interp(np.arange(x_max), BinCenters2, a_area2_pdf[1])
        fill_between(np.linspace(0,x_max,len(interp2)), interp2, interp1, facecolor='k', alpha=0.2)

        
        
        pl.plot_pdf(S_data, color='k', linewidth = line_width_fit, \
                                                           label='Before input')
    if regime == 'extrainput_start':
        
        a_area1_pdf = pl.pdf(a_area1, 10)
        a_area2_pdf = pl.pdf(a_area2, 10) 
        BinCenters1 = (a_area1_pdf[0][:-1]+a_area1_pdf[0][1:])/2.
        BinCenters2 = (a_area2_pdf[0][:-1]+a_area2_pdf[0][1:])/2.
        x_max = a_area1_pdf[0].max()
        interp1 = np.interp(np.arange(x_max), BinCenters1, a_area1_pdf[1])
        interp2 = np.interp(np.arange(x_max), BinCenters2, a_area2_pdf[1])
        fill_between(np.linspace(0,x_max,len(interp2)), interp2, interp1, facecolor='r', alpha=0.2)

        
        pl.plot_pdf(S_data, color='r', linewidth = line_width, \
                                                    label='Input onset')        
    if regime == 'extrainput_end':
        
        a_area1_pdf = pl.pdf(a_area1, 10)
        a_area2_pdf = pl.pdf(a_area2, 10) 
        BinCenters1 = (a_area1_pdf[0][:-1]+a_area1_pdf[0][1:])/2.
        BinCenters2 = (a_area2_pdf[0][:-1]+a_area2_pdf[0][1:])/2.
        x_max = a_area1_pdf[0].max()
        interp1 = np.interp(np.arange(x_max), BinCenters1, a_area1_pdf[1])
        interp2 = np.interp(np.arange(x_max), BinCenters2, a_area2_pdf[1])
        fill_between(np.linspace(0,x_max,len(interp2)), interp2, interp1, facecolor='cyan', alpha=0.2)

        
        pl.plot_pdf(S_data, color='cyan', linewidth = line_width_fit, \
                                                 label = 'Readaptation') 


subplot(121)
xscale('log'); yscale('log')
xlabel(r'$T$', fontsize=letter_size)
ylabel(r'$f(T)$', fontsize=letter_size) 
fig_6a.spines['right'].set_visible(False)
fig_6a.spines['top'].set_visible(False)
fig_6a.xaxis.set_ticks_position('bottom')
fig_6a.yaxis.set_ticks_position('left')

xlim([1, 300])
ylim([0.0001, 1])
xticks([1, 10, 100], ['$10^0$', '$10^{1}$', '$10^{2}$'])
yticks([1, 0.01, 0.0001], ['$10^0$', '$10^{-2}$', '$10^{-4}$'])
tick_params(labelsize=letter_size)                

              
subplot(122)
xscale('log'); yscale('log')	 
xlabel(r'$S$', fontsize=letter_size)
ylabel(r'$f(S)$', fontsize=letter_size)
fig_6b.spines['right'].set_visible(False)
fig_6b.spines['top'].set_visible(False)
fig_6b.xaxis.set_ticks_position('bottom')
fig_6b.yaxis.set_ticks_position('left')

xlim([1, 3000])
ylim([0.00001, 0.1])
xticks([1, 10, 100, 1000], \
     ['$10^0$', '$10^{1}$', '$10^{2}$', '$10^{3}$'])
yticks([0.1, 0.001, 0.00001],\
            ['$10^{-1}$', '$10^{-3}$', '$10^{-5}$'])
legend(loc=(0.5, 0.8), prop={'size':letter_size}, frameon=False)

########################################################################
fig_6a.annotate('A', xy=subplot_letter, xycoords='axes fraction', \
                fontsize=letter_size_panel ,  fontweight='bold', \
                horizontalalignment='right', verticalalignment='bottom') 
fig_6b.annotate('B', xy=subplot_letter, xycoords='axes fraction', \
                fontsize=letter_size_panel ,  fontweight='bold', \
                horizontalalignment='right', verticalalignment='bottom') 
                
gcf().subplots_adjust(bottom=0.17)
fig_6.subplots_adjust(wspace=.4) 
 
print 'Saving figures...',		
result_path = '../../../Avalanche_Results/'
result_name_png = 'Fig6_test.pdf'
savefig(os.path.join(result_path, result_name_png), format = 'pdf')
plt.show()
