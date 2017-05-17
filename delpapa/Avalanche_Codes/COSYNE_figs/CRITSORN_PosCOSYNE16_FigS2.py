####
# Script for the fifth paper figure
# Includes: Only a subset of neurons noise (A) and activity (B)  
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
THETA_P = 25           # here, half doesnt work (cause the activity is                                       

figure(1, figsize=(20,20))
fig_S2a = subplot(111)



########################################################################
# Fig. S2A,B: Fixed Spikes (power-laws and activity)

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
    T_data, S_data = analysis.avalanches(data_all, \
                                           'N', '200', Theta_percent=THETA_P)
                                           
    if p == '000':
        label_p = '$0% N^E$'
    elif p == '005':
        label_p = '$5% N^E$'
    elif p == '010':
        label_p = '$10% N^E$'
    pl.plot_pdf(S_data, linewidth=5.0, label=label_p) 
                                           
        
    xscale('log'); yscale('log')
    xlabel('S', fontsize=35); 
    ylabel('f(S)', fontsize=35);
    ylim([0.000001, 0.1])
    xlim([1, 3000])
    tick_params(axis='both', which='major', labelsize=35)
    legend(loc='best', prop={'size':35}, frameon=False)
        
        
########################################################################

print 'Saving figures...',		
result_path = '../Avalanche_Results/'
result_name_png = 'CRITSORN_PosCOSYNE16_FigS2.png'

figure(1)
savefig(os.path.join(result_path, result_name_png), format = 'png')
print 'done'

show()
