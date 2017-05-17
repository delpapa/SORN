####
# Script for the first paper figure
# Includes: 3 Network phases (A); avalanche definition (B) 
####


from pylab import * 

import tables
import os
from tempfile import TemporaryFile

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

figure(1, figsize=(15, 15))



## Avalanche definition
print 'Avalanche definition...'
########################################################################
# Fig. 1A: AVAL_DEF
fig_1a = subplot(111)

boundary = Theta*np.ones(plot_last_steps)
plot(activity, 'k', label='Network activity', linewidth=2.0)
plot(boundary, '--k', label='Activity threshold', linewidth=5.0)

fill_between(np.arange(plot_last_steps), activity, boundary, \
 alpha = 0.3, where=activity>=boundary, facecolor='r', interpolate=True)
 
# annotate stuff 
text(25, 45, r'avalanches', fontsize=35, color='k')
text(122, 10, r'T', fontsize=35, color='b') 
text(122, 20, r'S', fontsize=35, color='r') 
text(65, 1, r'100 time steps', fontsize=35, color='k')

#arrow(35, 44, -20, -16, head_width=3, head_length=5, fc='k', ec='k', length_includes_head=True, width=0.5)
#arrow(55, 44, 16, -9, head_width=3, head_length=5, fc='k', ec='k', length_includes_head=True, width=0.4)
#arrow(65, 44, 29, -9, head_width=3, head_length=5, fc='k', ec='k', length_includes_head=True, width=0.3)
                                                                        
                                                                                                                                  
annotate(s='', xy=(92,13), xytext=(157,13), \
                              arrowprops=dict(arrowstyle='<->', ec='b'))
annotate(s='', xy=(50,0), xytext=(150,0), \
                              arrowprops=dict(arrowstyle='-', ec='k'))

# axis stuff
xlim([0, plot_last_steps])
ylim([0, 50])

fig_1a.spines['right'].set_visible(False)
fig_1a.spines['top'].set_visible(False)
fig_1a.spines['bottom'].set_visible(False)
fig_1a.yaxis.set_ticks_position('left')
fig_1a.axes.get_xaxis().set_visible(False)

yticks(arange(0, 50, 20))

ylabel('a(t)', fontsize=35)
tick_params(axis='both', which='major', labelsize=35)

#legend stuff
legend(loc=(0.5, 0.85), prop={'size':35}, frameon=False)

########################################################################
# saving figures
print 'Saving figures...',		
result_path = '../Avalanche_Results/'
result_name_png = 'CRITSORN_PosCOSYNE16_Fig1.1.png'
   
figure(1)
savefig(os.path.join(result_path, result_name_png), format='png')
print 'done'
print '\n\n'

show()
