####
# Script for the 3rd paper figure (from PhengshengSORN - N200)
# A: activity distribution
# B: power-laws for different Theta
#
# Script also used for comparison between power-law exponents
# for different Theta for Table S5 (from Fig3B)
#
####

from pylab import *

import tables
import os
from tempfile import TemporaryFile
from matplotlib import gridspec

import data_analysis as analysis
import powerlaw as pl

exp_name = 'N200'
stable_steps = 3e6
number_of_files = 50
Theta_range = np.array([5, 10, 15, 20, 25]) # values of Theta

### figure parameters
width  =  8
height = 3
fig = figure(1, figsize=(width, height))
gs = gridspec.GridSpec(1, 2)
letter_size = 10
letter_size_panel = 12
line_width = 1.5
line_width_fit = 2.0
subplot_letter = (-0.25, 1)

### color parameters
c_size = '#B22400'
c_duration = '#006BB2'
c_rawdata = 'gray'
c_expcut = 'k'
c_stable = '#2E4172'
c_notstable = '#7887AB'
########################################################################
# Fig. 3A: activity distribution

print 'Fig. 3A...'
fig_3a = subplot(121)

data_all = np.zeros((number_of_files, stable_steps))
for result_file in range(number_of_files):

    exper_path = ''
    h5 = tables.openFile(os.path.join(exper_path,'result.h5'),'r')
    data = h5.root
    data_all[result_file] = \
                  np.around(data.activity[0][-stable_steps:]*data.c.N_e)
    h5.close()

act_density = zeros((number_of_files, data_all.max()+1))
for data_file in xrange(number_of_files):

    print 'Activity from file ' + str(data_file+1)

    for i in xrange(int(stable_steps)):
        act_density[data_file, data_all[data_file, i]] += 1
    act_density[data_file, :] /= act_density[data_file, :].sum()

# activity distribution and std
act_density_mean = act_density.mean(0)
act_density_std = act_density.std(0)

plot(act_density_mean, c_stable, linewidth=line_width)


text(5, 0.062, r'$5\%$', fontsize=letter_size, color='k')
text(14, 0.062, r'$25\%$', fontsize=letter_size, color='k')
text(6, 0.072, r'$\langle a(t) \rangle_t / 2$', \
                                        fontsize=letter_size, color='k')

# axis stuff
tick_params(axis='both', which='major', labelsize=letter_size)
ylabel(r'$p(a(x))$', fontsize=letter_size)
xlabel(r'$a(x)$ [# neurons]', fontsize=letter_size)
xlim(0, act_density_mean.size)
ylim(0, 0.07)

xticks([0, 20, 40, 60], \
       ['$0$', '$20$', '$40$', '$60}$'])
yticks([0, 0.03, 0.06], \
       ['$0.00$', '$0.03$', '$0.06$'])

fig_3a.spines['right'].set_visible(False)
fig_3a.spines['top'].set_visible(False)
fig_3a.tick_params(axis=u'both', which=u'both',length=0)
########################################################################

########################################################################
# Fig. 3B: avalanche for different Theta

print 'Fig. 3B...'
fig_3b = subplot(122)

print data_all.mean()

for Theta_percent in Theta_range:
    print 'Theta = '  + str(Theta_percent) + '%'

    # Theta percentile goes here
    T_data, S_data = analysis.avalanches(data_all, 'N', '200', \
					                      Theta_percent = Theta_percent)

    pl.plot_pdf(S_data, label =  str(Theta_percent)+ r'%', \
                                                   linewidth=line_width)

    ####################################################################
    # alpha and tau exponents calculation for Table S5
    # takes time - comment out if not needed
    T_fit = pl.Fit(T_data, xmin=6, xmax=60, discrete=True)
    T_alpha = T_fit.alpha
    T_sigma = T_fit.sigma
    S_fit = pl.Fit(S_data, xmin=10, xmax=1500, discrete=True)
    S_alpha = S_fit.alpha
    S_sigma = S_fit.sigma
    print 'alpha = ', T_alpha, 'sigma = ', T_sigma
    print 'Loglikelyhood ratio (power-law/exp) =', \
               T_fit.distribution_compare('power_law','exponential', \
                        normalized_ratio=True)
    print 'tau = ', S_alpha, 'sigma = ', S_sigma
    print 'Loglikelyhood ratio (power-law/exp) =', \
               S_fit.distribution_compare('power_law','exponential', \
                        normalized_ratio=True)
    ####################################################################

    if Theta_percent == Theta_range.min():
        theta_low = percentile(data_all, Theta_percent)
        fig_3a.plot((theta_low, theta_low), (0, 0.06), 'k--')
    if Theta_percent == Theta_range.max():
        theta_high = percentile(data_all, Theta_percent)
        fig_3a.plot((theta_high, theta_high), (0, 0.06), 'k--')
    if Theta_percent == 10:
        theta_aver = percentile(data_all, data_all.mean()/2)
        fig_3a.plot((theta_aver, theta_aver), (0, 0.07), 'k--', \
                                               linewidth=line_width_fit)

### Fill the shaded area
range_to_fill = arange(0, len(act_density_mean), 1)
wh_to_fill = (range_to_fill <= theta_high) & (range_to_fill >=theta_low)
boundary = zeros(len(act_density_mean))
fig_3a.fill_between(range_to_fill, act_density_mean, boundary,\
                      where = wh_to_fill, facecolor=c_stable, alpha=0.5)

### axis stuff
xscale('log'); yscale('log')
xlabel(r'$S$', fontsize=letter_size)
ylabel(r'$f(S)$', fontsize=letter_size)

fig_3b.spines['right'].set_visible(False)
fig_3b.spines['top'].set_visible(False)
fig_3b.xaxis.set_ticks_position('bottom')
fig_3b.yaxis.set_ticks_position('left')
tick_params(labelsize=letter_size)

# ticks name
xlim([1, 3000])
ylim([0.00001, 0.1])
xticks([1, 10, 100, 1000], \
     ['$10^0$', '$10^{1}$', '$10^{2}$', '$10^{3}$'])
yticks([0.1, 0.001, 0.00001],\
            ['$10^{-1}$', '$10^{-3}$', '$10^{-5}$'])

# legend stuff
legend(loc=(0.5, 0.55), prop={'size':letter_size}, \
                            title=r'Activity percentile', frameon=False)
fig_3b.get_legend().get_title().set_fontsize(letter_size)


########################################################################
fig_3a.annotate('A', xy=subplot_letter, xycoords='axes fraction', \
                fontsize=letter_size_panel ,  fontweight='bold', \
                horizontalalignment='right', verticalalignment='bottom')
fig_3b.annotate('B', xy=subplot_letter, xycoords='axes fraction', \
                fontsize=letter_size_panel ,  fontweight='bold', \
                horizontalalignment='right', verticalalignment='bottom')
gcf().subplots_adjust(bottom=0.17)
fig.subplots_adjust(wspace=.4)

print 'Saving figures...',
result_path = '../../plot'
result_name_png = 'Fig3.pdf'
savefig(os.path.join(result_path, result_name_png), format='pdf')
print 'done\n\n'
