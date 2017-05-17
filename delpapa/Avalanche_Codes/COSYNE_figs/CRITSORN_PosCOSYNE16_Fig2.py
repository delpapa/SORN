####
# Script for the second paper figure
# Includes: N200 duration (A); N200 size (B)
#           Different network sizes duration (C) and size (D)
#           Ratio between exponents (E)
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

figure(1, figsize=(25,13))
########################################################################
# Fig. 1A and 1B SORN size and duration with exponents (A and B)

exp_name = 'N200'
number_of_files = 50
stable_steps = 3e6 # steps to use: after transient (2e6)
THETA = 'half'

#count files
print 'Loading experiment files...'
exper_path = '../Avalanche_Experiments/Pengsheng_SORN/' + exp_name

#load the data in 'data_all'	
data_all = zeros((number_of_files, stable_steps))    
for result_file in range(number_of_files):
        result_path = exper_path +'/'+str(result_file+1)+'/common/'
        h5 = tables.openFile(os.path.join(result_path,'result.h5'),'r')
        data = h5.root
        data_all[result_file] = \
                     around(data.activity[0][-stable_steps:]*data.c.N_e) 
        h5.close()

# calculates avalanches    
T_data, S_data = analysis.avalanches(data_all, \
                            exp_name[0], exp_name[1:], Threshold=THETA)
                            
########################################
### duration stuff
fig_1a = subplot(122)
print 'Fig. 1A...'

# raw data
T_x, inverse = unique(T_data, return_inverse=True)
y_freq = bincount(inverse)
T_y = y_freq / float(y_freq.sum()) # normalization
plot(T_x, T_y, 'rx')

# power law fit
T_fit = pl.Fit(T_data, xmin=6, xmax=80, discrete=True)
T_alpha = T_fit.alpha
T_sigma = T_fit.sigma
T_xmin = T_fit.xmin
#pl.plot_pdf(T_data, color='b', linewidth=4.0)

T_fit.power_law.plot_pdf(color='k', \
            label = r'power law: $ \alpha = $' + str(round(T_alpha, 2)), linewidth=5.0)
            
# exp cutoff calculation
T_fit = pl.Fit(T_data, xmin=6, discrete=True)#, xmin_distance='Asquare')
T_trunc_alpha = T_fit.truncated_power_law.parameter1
T_trunc_beta = T_fit.truncated_power_law.parameter2

T_fit.truncated_power_law.plot_pdf(color='b', \
     label = r'exp. cut-off: $ \alpha =$' + str(round(T_trunc_alpha, 2)), linewidth=5.0) 

# axis stuff
xscale('log'); yscale('log')
xlabel('T', fontsize=35)
ylabel('f(T)', fontsize=35)
xlim(1, 1000)
ylim([0.00001, 1])
tick_params(axis='both', which='major', labelsize=35)

# legend stuff
legend(loc='best', prop={'size':35}, frameon=False)

#########################################
### size stuff
fig_1b = subplot(121)
print 'Fig. 1B...'

# raw data
S_x, inverse = unique(S_data, return_inverse=True)
y_freq = bincount(inverse)
S_y = y_freq / float(y_freq.sum())
plot(S_x, S_y, 'rx')

# power law fit
S_fit = pl.Fit(S_data, xmin=10, xmax=1200, discrete=True)
S_alpha = S_fit.alpha
S_sigma = S_fit.sigma
S_xmin = S_fit.xmin
#pl.plot_pdf(S_data, color='b', linewidth=4.0) 

S_fit.power_law.plot_pdf(color='k', \
            label = r'power law: $ \tau = $' + str(round(S_alpha, 2)), linewidth=5.0)

                             
# exp cutoff calculation
S_fit = pl.Fit(S_data, xmin=10, discrete=True)
S_trunc_alpha = S_fit.truncated_power_law.parameter1
S_trunc_beta = S_fit.truncated_power_law.parameter2


S_fit.truncated_power_law.plot_pdf(color='b', \
    label = r'exp. cut-off: $ \tau =$' + str(round(S_trunc_alpha, 2)), linewidth=5.0)  


# axis stuff
xscale('log'); yscale('log')
xlabel('S', fontsize=35)
ylabel('f(S)', fontsize=35)
xlim([1, 10000])
ylim([0.000001, 0.1])
tick_params(axis='both', which='major', labelsize=35)

# legend stuff
legend(loc='best', prop={'size':35}, frameon=False)


# saving figures
print 'Saving figures...',		
result_path = '../Avalanche_Results/'
result_name_png = 'CRITSORN_PosCOSYNE16_Fig2.png'
   
figure(1)
savefig(os.path.join(result_path, result_name_png), format='png')
print 'done'
print '\n\n'
