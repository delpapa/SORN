####
# Script for the fifth paper figure
# Includes: Gaussian noise (A) and activity (B)
#           Random neuron noise (C) and activity (D)
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
stable_steps = 3e6
number_of_files = 50
THETA = 'half'

figure(1, figsize=(25,12))
fig_5a = subplot(121)
fig_5b = subplot(122)


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
        exper_path =  '../Avalanche_Experiments/PhaseTransition_Noise/'\
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
    subplot(121)
    pl.plot_pdf(S_data, linewidth=5.0)
    
    if sigma == '0.05':
        # power law fit
        S_fit = pl.Fit(S_data, xmin=10, xmax=1500, discrete=True)
        S_alpha = S_fit.alpha
        S_xmin = S_fit.xmin
        S_fit.power_law.plot_pdf(color='k', \
            label = r'power law: $ \tau = $' + str(round(S_alpha, 2)), \
                                                          linewidth=5.0)
    legend(loc='best', prop={'size':35}, frameon=False)
    
        
    xscale('log'); yscale('log')
    xlabel('S', fontsize=35); 
    ylabel('f(S)', fontsize=35);
    ylim([10e-10, 10e-1])
    xlim(1, 100000)
    ylim([0.00000001, 1])
    tick_params(axis='both', which='major', labelsize=35)
        
        
    subplot(122) 
    plot(act_density_mean, linewidth=5.0, label=r'$\sigma^2 = $'+sigma) 
        
    ylabel('Activity distribution', fontsize=35)
    xlabel('# Active neurons', fontsize=35)
    xlim([0, 70])
    yticks(arange(0, 0.05, 0.11))
    legend(loc='best', prop={'size':35}, frameon=False)
    tick_params(axis='both', which='major', labelsize=35)

########################################################################

print 'Saving figures...',		
result_path = '../Avalanche_Results/'
result_name_png = 'CRITSORN_PosCOSYNE16_Fig5.png'

figure(1)
savefig(os.path.join(result_path, result_name_png), format = 'png')
print 'done'

show()
