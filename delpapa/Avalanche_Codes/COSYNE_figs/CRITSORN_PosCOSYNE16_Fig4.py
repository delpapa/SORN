####
# Script for the fourth paper figure
# Includes: Frozen plasticity comparison (A)
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
number_of_files = 50    
THETA = 'half'                                                

figure(1, figsize=(25, 13))
fig_4a = subplot(122)
fig_4b = subplot(121)


########################################################################
# Fig. 4A,B: Compare distributions of AllFrozen and normal SORN

for experiment_name in ['AllFrozen']:

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

        subplot(121)
        if regime == 'non-frozen':
            pl.plot_pdf(a_dur, color='k', linewidth=5.0)
        if regime == 'frozen':
            pl.plot_pdf(a_dur, linewidth=5.0)
            
        subplot(122)
        if regime == 'non-frozen':
            pl.plot_pdf(a_area, color='k', linewidth=5.0, label ='SORN')
        if regime == 'frozen':             
            pl.plot_pdf(a_area, linewidth=5.0, label='Frozen (All)')


transient_steps = 0 # special way to plot an only initialized SORN
print 'Freezing plasticity in the beginning...'
result_path ='../Avalanche_Results/Frozen_Plasticity/' + \
                                                      'AllFrozen_Step0/'                                                   
for regime in ['frozen']:
        
    print '\n', regime
    data_all = np.zeros((number_of_files, frozen_steps))    
        
    for result_file in range(number_of_files):
        exper = 'result.h5'
        exper_path =  '../Avalanche_Experiments/Frozen_Plasticity/'\
                  + 'AllFrozen_Step0/' + str(result_file+1) + '/common/'
        h5 = tables.openFile(os.path.join(exper_path,exper),'r')
        data = h5.root
        
        data_all[result_file] = np.around(data.activity[0] \
                             [frozen_steps+transient_steps:]*data.c.N_e)            
        
        h5.close()
            
    a_dur, a_area = analysis.avalanches(data_all, 'N', '200', \
                                                       Threshold=THETA)
        
    subplot(121)
    pl.plot_pdf(a_dur, color='r', linewidth=5.0)
            
    subplot(122)
    pl.plot_pdf(a_area, color='r', linewidth=5.0, \
                                                 label='Random network')
                                                 
# axis stuff    
subplot(121)
xscale('log'); yscale('log')
xlim(0, 1e4)
xlabel('T', fontsize=35)
ylabel('f(T)', fontsize=35)
xlim(1, 10000)
ylim([0.00000001, 1])
tick_params(axis='both', which='major', labelsize=35)

subplot(122)
xscale('log'); yscale('log')
xlim(0, 1e5)	 
xlabel('S', fontsize=35) 
ylabel('f(S)', fontsize=35)
xlim(1, 100000)
ylim([0.000000001, 0.1])
legend(loc='best', prop={'size':35}, frameon=False)
tick_params(axis='both', which='major', labelsize=35)
########################################################################    
    
print 'Saving figures...',		
result_path = '../Avalanche_Results/'
result_name_png = 'CRITSORN_PosCOSYNE16_Fig4.png'
figure(1)
savefig(os.path.join(result_path, result_name_png), format = 'png')
print 'done'

show()
