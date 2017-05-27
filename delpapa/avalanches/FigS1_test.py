####
# Script for the S1 paper figure (from FrozenPlasticity)
# A: avalanches for S for full area
# B: avalanches for different bin sizes
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

number_of_files = 50
THETA = 'half'                                                

width  =  9
height = 5 # fix this to the golden ratio
fig_S1 = figure(1, figsize=(width, height))
fig_S1a = subplot(121)
fig_S1b = subplot(122)

letter_size = 14
letter_size_panel = 14
line_width = 1.5
line_width_fit = 2.0
subplot_letter = (-0.2, 1)

########################################################################
# Fig. B: bins sizes
fig_S1b = subplot(122)
stable_steps = 3e6 # steps to use: after transient (2e6)
THETA = 'half'

#count files
print 'Loading experiment files...'
exper_path = '../Avalanche_Experiments/Pengsheng_SORN/N200'

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
for binsize in [1, 2, 10, 20]:
    
    T_data, S_data = analysis.avalanches(data_all, \
             'N', '200', Threshold=THETA)

    # power law fit
    S_fit = pl.Fit(S_data, xmin=10, xmax=1500, discrete=True)
    S_alpha = S_fit.alpha
    S_sigma = S_fit.sigma
    S_xmin = S_fit.xmin
    pl.plot_pdf(S_data, binsize=binsize, linewidth=line_width, label = r'$b_s = $'+str(1./binsize))

### axis stuff
xscale('log'); yscale('log')
xlabel(r'$S$', fontsize=letter_size)
ylabel(r'$f(S)$', fontsize=letter_size)
tick_params(axis='both', which='major', labelsize=letter_size)

fig_S1b.spines['right'].set_visible(False)
fig_S1b.spines['top'].set_visible(False)
fig_S1b.tick_params(axis=u'both', which=u'both',length=0)

# ticks name
xlim([1, 3000])
ylim([0.00001, 0.1])
xticks([1, 10, 100, 1000], \
     ['$10^0$', '$10^{1}$', '$10^{2}$', '$10^{3}$'])
yticks([0.1, 0.001, 0.00001],\
            ['$10^{-1}$', '$10^{-3}$', '$10^{-5}$'])

# legend stuff
legend(loc=(0.45, 0.65), prop={'size':letter_size}, frameon=False)
fig_S1b.get_legend().get_title().set_fontsize(letter_size)

########################################################################
tight_layout()
fig_S1a.annotate('A', xy=subplot_letter, xycoords='axes fraction', \
                fontsize=letter_size_panel,  fontweight='bold', \
                horizontalalignment='right', verticalalignment='bottom') 
fig_S1b.annotate('B', xy=subplot_letter, xycoords='axes fraction', \
                fontsize=letter_size_panel,  fontweight='bold', \
                horizontalalignment='right', verticalalignment='bottom') 

print 'Saving figures...',		
result_path = '../Avalanche_Results/'
result_name_png = 'FigS1_test.pdf'
savefig(os.path.join(result_path, result_name_png), format = 'pdf')
