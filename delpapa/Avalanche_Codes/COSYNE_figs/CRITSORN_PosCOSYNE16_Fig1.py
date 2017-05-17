####
# Script for the first COSYNE16 figure
# Includes: Self-organization SORN phases 
####


from pylab import * 

import tables
import os
from tempfile import TemporaryFile

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

import matplotlib.ticker as mtick

variable = 'N'; value = '200'
exp_file = 19 # up to 50
plot_last_steps = 235 

exper = 'result.h5'
exper_path = '../Avalanche_Experiments/Pengsheng_SORN/' \
                   + variable + value + '/' + str(exp_file) + '/common/'
h5 = tables.openFile(os.path.join(exper_path,exper),'r')
data = h5.root	
activity = data.activity[0][-plot_last_steps-50:-50]*float(value)
connec_frac = data.ConnectionFraction[0]
Theta = activity.mean() / 1.3
h5.close()	

figure(1, figsize=(35, 12))



########################################################################
# Fig. 1: CONNEC_FRAC
fig_1 = subplot(111)
plot(connec_frac[:2e6]*100, 'k', linewidth=5)
plot(arange(2e6, 3e6), connec_frac[2e6:3e6]*100, 'r', linewidth=5)

#annotate stuff
text(2e6, 2, r'decay', fontsize=35, color='k')
text(0.6e6, 10, r'growth', fontsize=35, color='k') 
text(25e5, 10, r'stable', fontsize=35, color='r')  

fig_1.spines['right'].set_visible(False)
fig_1.spines['top'].set_visible(False)
fig_1.yaxis.set_ticks_position('left')
fig_1.xaxis.set_ticks_position('bottom')

xticks(arange(0, 3.1e6, 1e6))
yticks(arange(0, 17, 5))

fmt = '%2.0f%%'
yticks = mtick.FormatStrFormatter(fmt)
fig_1.yaxis.set_major_formatter(yticks)

formatter = matplotlib.ticker.ScalarFormatter()
formatter.set_powerlimits((-3,4))
fig_1.xaxis.set_major_formatter(formatter)

xlabel(r'$10^6$ time steps', fontsize=35)
ylabel('Active Connections', fontsize=35)
tick_params(axis='both', which='major', labelsize=35)
########################################################################

# saving figure
print 'Saving figures...',		
result_path = '../Avalanche_Results/'
result_name_png = 'CRITSORN_PosCOSYNE16_Fig1.png'
   
figure(1)
savefig(os.path.join(result_path, result_name_png), format='png')
print 'done'
print '\n\n'

show()
