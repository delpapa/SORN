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
number_of_files = 50
THETA = 'half'
raster_steps = 800

### figure parameters
width  =  6
height = 4.2
fig = figure(1, figsize=(width, height))
    
fig_5a = subplot(3, 1, 1)
fig_5b = subplot(3, 1, 2)
fig_5c = subplot(3, 1, 3)

letter_size = 11
letter_size_panel = 13
line_width = 1.5
line_width_fit = 2.0
subplot_letter = (-0.18, 1.05)
subplot_letter_long = (-0.08, 0.85)

######################################################################
# Fig. 5E, F, G: Raster plots

experiment_name = 'Gaussian'
for sigma in ['0.005', '0.05', '5']:
        
    for data_file in xrange(1):
        exper = 'result.h5'
        exper_path =  '../../../Avalanche_Experiments/PhaseTransition_Noise/'\
         + experiment_name + '/sigma' + sigma + '/' + str(data_file+1)\
                                                            + '/common/'
        h5 = tables.openFile(os.path.join(exper_path,exper),'r')
        data = h5.root
    
        ### Save raster plots
        if data_file == 0: # 0 is a random file
            if sigma == '0.005':
                raster_low = data.Spikes[0][:, -raster_steps:]
            if sigma == '0.05':
                raster_inter = data.Spikes[0][:, -raster_steps:]
            if sigma == '5':
                raster_high = data.Spikes[0][:, -raster_steps:]                            
        h5.close()

subplot(311)
for (i,sp) in enumerate(raster_low):
    s_train = where(sp == 1)[0]
    if s_train != []:
        vlines(s_train, i + 0.5, i + 1.5)
        hold('on')
ylabel('Low', fontsize=letter_size)
ylim([0, 200])
tick_params(axis='both', which='major', labelsize=letter_size) 
xticks([])
yticks([])

subplot(312)
for (i,sp) in enumerate(raster_inter):
    s_train = where(sp == 1)[0]
    if s_train != []:
        vlines(s_train, i + 0.5, i + 1.5)
        hold('on')
ylabel('Interm.', fontsize=letter_size)
ylim([0, 200])
tick_params(axis='both', which='major', labelsize=letter_size)
xticks([]) 
yticks([])  

subplot(313)
for (i,sp) in enumerate(raster_high):
    s_train = where(sp == 1)[0]
    if s_train != []:
        vlines(s_train, i + 0.5, i + 1.5)
        hold('on')
xlabel('Time step', fontsize=letter_size)
ylabel('High', fontsize=letter_size)
ylim([0, 200]) 
tick_params(axis='both', which='major', labelsize=letter_size) 
yticks([]) 
xticks([0, 400, 800], \
       ['$0$', '$400$', '$800$'])


########################################################################
fig_5a.annotate('E', xy=subplot_letter_long, xycoords='axes fraction', \
                fontsize=letter_size_panel ,  fontweight='bold', \
                horizontalalignment='right', verticalalignment='bottom') 
fig_5b.annotate('F', xy=subplot_letter_long, xycoords='axes fraction', \
                fontsize=letter_size_panel ,  fontweight='bold', \
                horizontalalignment='right', verticalalignment='bottom')
fig_5c.annotate('G', xy=subplot_letter_long, xycoords='axes fraction', \
                fontsize=letter_size_panel ,  fontweight='bold', \
                horizontalalignment='right', verticalalignment='bottom')                 
fig.subplots_adjust(wspace=.25)  
fig.subplots_adjust(hspace=.45)

print 'Saving figures...',		
result_path = '../../../Avalanche_Results/'
result_name_png = 'Fig5_2.eps'
savefig(os.path.join(result_path, result_name_png), format='eps')
print 'done\n\n'
