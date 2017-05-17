####
# Script for the sixth paper figure
# Includes: noise-ABCD-noise
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
sequence_steps = 2e6
number_of_files = 20
experiment_folder = 'Noise_Random(0.1)NoPlasticity_Random(0.1)Plasticity'
possible_regimes = ['noise1', 'sequence', 'noise2']


figure(1, figsize=(25, 13))

########################################################################
# Fig. 1A,B: Gaussian Noise + ABCD + Gaussian Noise

for regime in possible_regimes:
    
    print '\n', regime, '...'
    data_all = zeros((number_of_files, section_steps))
    if regime == 'sequence': 
        data_all = zeros((number_of_files, sequence_steps))
    #####
    act_density = zeros((number_of_files, 200)) 
    #####
    for result_file in xrange(number_of_files):
    
        exper = 'result.h5'
        exper_path =  '../Avalanche_Experiments/Extra_Input/' + \
                                                   experiment_folder + \
                                   '/' + str(result_file+1) + '/common/'
        h5 = tables.openFile(os.path.join(exper_path,exper),'r')
        data = h5.root
        
        if regime == 'noise1':
            data_all[result_file] = np.around(data.activity[0] \
                             [section_steps:2*section_steps]*data.c.N_e)
        if regime == 'sequence':
            
            data_all[result_file] = np.around(data.activity[0] \
                           [2*section_steps:2*section_steps+sequence_steps]*data.c.N_e)                                                                                              
        if regime == 'noise2':
            data_all[result_file] = np.around(data.activity[0] \
                           [2*section_steps+sequence_steps:3*section_steps+sequence_steps]*data.c.N_e)  
        
        #### 
        for i in xrange(int(section_steps)):
            act_density[result_file, data_all[result_file, i]] += 1
        act_density[result_file, :] /= act_density[result_file, :].sum() 
        ####      
                                             
        h5.close()

    print 'Mean =', data_all.mean()
    print 'Std =', data_all.std()
    
    #################################################
    figure(2)
    act_density = zeros((number_of_files, data_all.max()+1))    
    for data_file in xrange(number_of_files):
        
        for i in xrange(int(section_steps)):
            act_density[data_file, data_all[data_file, i]] += 1
        act_density[data_file, :] /= act_density[data_file, :].sum()
    
    # activity distribution   
    act_density_mean = act_density.mean(0)
    if regime == 'noise1':
        plot(act_density_mean, 'b')
    if regime == 'sequence':
        plot(act_density_mean, 'r')
    if regime == 'noise1':
        plot(act_density_mean, 'g')    
    figure(1)
    #################################################################
    
    if regime == 'noise1':
        THETA = int(data_all.mean()/2)
    
    if regime in ['noise1','noise2']:    
        T_data, S_data = analysis.avalanches(data_all, 'N', '200', \
                                                        Threshold=10)
    else:
        T_data, S_data = analysis.avalanches(data_all, 'N', '200', \
                                                     #Theta_percent=21)
                                                     Threshold=22)
    
    #~ if regime == 'noise1':
        #~ 
        #~ T_fit = pl.Fit(T_data, xmin=6, xmax=80, discrete=True)
        #~ T_alpha = T_fit.alpha
        #~ subplot(122)
        #~ T_fit.power_law.plot_pdf(color='k', linewidth=5.0)
#~ 
        #~ S_fit = pl.Fit(S_data, xmin=10, xmax=1200, discrete=True)
        #~ S_alpha = S_fit.alpha
        #~ subplot(121)
        #~ S_fit.power_law.plot_pdf(color='k', \
          #~ label = r'$\alpha = $' + str(round(T_alpha, 2)) + \
                  #~ r',$\tau = $' + str(round(S_alpha, 2)) , linewidth=5.0)

    ######
    subplot(122)
    if regime == 'noise1':
        pl.plot_pdf(T_data, color='b', linewidth = 5.0, label='Normal SORN')
    if regime == 'sequence':
        pl.plot_pdf(T_data, color='r', linewidth = 5.0, \
                                                    label='Extra input (no plasticity)')         
    if regime == 'noise2':
        pl.plot_pdf(T_data, color='g', linewidth = 5.0, label = 'Extra input (with plasticity)')

    subplot(121)
    if regime == 'noise1':
        pl.plot_pdf(S_data, color='b', linewidth = 5.0)
    if regime == 'sequence':
        pl.plot_pdf(S_data, color='r', linewidth = 5.0)        
    if regime == 'noise2':
        pl.plot_pdf(S_data, color='g', linewidth = 5.0) 



subplot(121)
xscale('log'); yscale('log')	 
xlabel('S', fontsize=30) 
ylabel('f(S)', fontsize=30)
xlim(1, 10000)
ylim([0.000001, 0.1])
tick_params(axis='both', which='major', labelsize=35)

subplot(122)
xscale('log'); yscale('log')
xlabel('T', fontsize=30); 
ylabel('f(T)', fontsize=30) 
xlim(1, 1000)
ylim([0.00001, 1])
tick_params(axis='both', which='major', labelsize=35)
legend(loc=(0.19, 0.8), prop={'size':35}, frameon=False) 

########################################################################
	
print 'Saving figures...',		
result_path = '../Avalanche_Results/'
result_name_png = 'CRITSORN_PosCOSYNE16_Fig6_new.png'

figure(1)
savefig(os.path.join(result_path, result_name_png), format = 'png')
print 'done'

show()






