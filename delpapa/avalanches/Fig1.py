################################################################################
# Script for the 1st paper figure                                              #
# A: 3 Network phases;                                                         #
# B: activity thresholding example                                             #
################################################################################

from pylab import *
import tables
import os
from tempfile import TemporaryFile

import matplotlib.ticker as mtick
import matplotlib.patches as patches
from matplotlib import gridspec


plot_last_steps = 150
exper = 'result.h5'

# this should be the path to the experiments folder
# CHANGE ACCORDINGLY!
# path for the experiment
exper_path = ''
# path for the results
result_path = '../../plots/'
h5 = tables.openFile(os.path.join(exper_path,exper),'r')
data = h5.root
activity = data.activity[0][-plot_last_steps-100:-100]*float(value)
connec_frac = data.ConnectionFraction[0]
Theta = 10
h5.close()

### figure parameters
width  =  10
height = 3 # width / 1.618 # fix this to the golden ratio
fig = figure(1, figsize=(width, height))
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 0.7])
letter_size = 13
line_width = 1.0

c_size = '#B22400'
c_duration = '#006BB2'
c_stable = '#2E4172'
c_notstable = '#7887AB'

################################################################################
# Fig. 1A: CONNEC_FRAC
fig_1a = subplot(gs[0])

### plot stuff
plot(connec_frac[:2e6]*100, c_notstable, linewidth=line_width)
plot(arange(2e6, 4e6), connec_frac[2e6:4e6]*100, c_stable, linewidth=line_width)

### annotate stuff
text(2e5, 1, r'decay', fontsize=letter_size, color=c_notstable)
text(8e5, 10, r'growth', fontsize=letter_size, color=c_notstable)
text(26e5, 13.5, r'stable', fontsize=letter_size, color=c_stable)

# axis stuff
xlim([0, 4e6])
ylim([0.1, 0.20])

fig_1a.spines['right'].set_visible(False)
fig_1a.spines['top'].set_visible(False)
fig_1a.spines['left'].set_visible(False)
fig_1a.spines['bottom'].set_visible(False)
fig_1a.yaxis.set_ticks_position('left')
fig_1a.xaxis.set_ticks_position('bottom')
fig_1a.tick_params(axis=u'both', which=u'both',length=0)
fig_1a.grid()

xticks(arange(0, 4.1e6, 1e6), ['0', '1', '2', '3', '4'])
yticks([5, 10, 15], ['5%', '10%', '15%'])

xlabel(r'$10^6$ time steps', fontsize=letter_size)
ylabel(r'Active Connections', fontsize=letter_size)
tick_params(axis='both', which='major', labelsize=letter_size)
################################################################################

################################################################################
### Fig. 1B: AVAL_DEF
fig_1b = subplot(gs[1])

### plot stuff
boundary = Theta*np.ones(plot_last_steps)
plot(activity, 'k', label='network activity', linewidth=line_width)
plot(boundary, '--k', label='$\theta$', linewidth=line_width)
fill_between(np.arange(plot_last_steps), activity, boundary, \
      alpha = 0.5, where=activity>=boundary, facecolor=c_size, interpolate=True)

### annotate stuff
text(20, 45, r'avalanches', fontsize=letter_size, color='k')
text(70, 4, r'duration', fontsize=letter_size, color=c_duration)
text(80, 12, r'size', fontsize=letter_size, color=c_size)
text(62, -4, r'100 time steps', fontsize=letter_size, color='k')

plot((58, 122), (8, 8) , c_duration, linewidth=2.0)
plot((50, 150), (0, 0) , 'k', linewidth=2.5)

### arrow stuff
arrow1 = patches.FancyArrowPatch((35,44), (12,29), arrowstyle='-|>', \
                                                fc='k', lw=1, mutation_scale=10)
fig_1b.add_patch(arrow1)
arrow2 = patches.FancyArrowPatch((55,44), (40,35), arrowstyle='-|>', \
                                                fc='k', lw=1, mutation_scale=10)
fig_1b.add_patch(arrow2)
arrow3 = patches.FancyArrowPatch((65,44), (75,38), arrowstyle='-|>', \
                                                fc='k', lw=1, mutation_scale=10)
fig_1b.add_patch(arrow3)

arrow4 = patches.FancyArrowPatch((60,44), (54.5,15), arrowstyle='-|>', \
                                                fc='k', lw=1, mutation_scale=10)
fig_1b.add_patch(arrow4)

### axis stuff
xlim([0, plot_last_steps])
ylim([0, 50])

fig_1b.spines['right'].set_visible(False)
fig_1b.spines['top'].set_visible(False)
fig_1b.spines['bottom'].set_visible(False)
fig_1b.yaxis.set_ticks_position('left')
fig_1b.axes.get_xaxis().set_visible(False)

yticks([0, Theta, 20, 40], ['0', r'$\theta$', '20', '40'])

ylabel(r'$a(t)$' + r' [# neurons]', fontsize=letter_size)
tick_params(axis='both', which='major', labelsize=letter_size)
################################################################################

### Panel label stuff
fig.text(0.01, 0.9, "A", weight="bold", \
                         horizontalalignment='left', verticalalignment='center')
fig.text(0.55, 0.9, "B", weight="bold", \
                         horizontalalignment='left', verticalalignment='center')
gcf().subplots_adjust(bottom=0.17)
fig.subplots_adjust(wspace=.4)

print 'Saving figures...',
result_name_png = 'Fig1.pdf'
savefig(os.path.join(result_path, result_name_png), format='pdf')
print 'done\n\n'
