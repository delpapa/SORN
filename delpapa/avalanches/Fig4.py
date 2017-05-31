####
# Script for the 4th paper figure (from FrozenPlasticity)
# A, B: AllFrozen and AllFrozen_Step0 (random network)
# C, D: IPFrozen and AllbutIPFrozen
####

from pylab import *

import tables
import os
from tempfile import TemporaryFile

import data_analysis as analysis

# work around to run powerlaw package [Alstott et al. 2014]
import powerlaw as pl

transient_steps = 2e6
frozen_steps = 2e6
number_of_files = 50
THETA = 'half'

width  =  6
height = 6
fig = figure(1, figsize=(width, height))
fig_4a = subplot(221)
fig_4b = subplot(222)
fig_4c = subplot(223)
fig_4d = subplot(224)

letter_size = 10
letter_size_panel = 12
line_width = 1.25
line_width_fit = 2.0
subplot_letter = (-0.2, 1)

### color parameters
c_size = '#B22400'
c_duration = '#006BB2'
c_rawdata = 'gray'
c_expcut = 'k'
c_stable = '#2E4172'
c_notstable = '#7887AB'

########################################################################
# Fig. 4A,B: Compare distributions of AllFrozen and normal SORN

for experiment_name in ['AllFrozen']:

    print experiment_name
    exper_path =''

    for regime in ['frozen', 'non-frozen']:

        print '\n', regime
        data_all = np.zeros((number_of_files, frozen_steps))

        for result_file in range(number_of_files):
            exper = 'result.h5'
            exper_path =  '../../../Avalanche_Experiments/Frozen_Plasticity/'\
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

        a_dur1, a_area1 = analysis.avalanches(data_all, 'N', '200',\
                                                      Theta_percent = 5)
        a_dur2, a_area2 = analysis.avalanches(data_all, 'N', '200',\
                                                     Theta_percent = 25)



        subplot(221)
        if regime == 'non-frozen':

            pl.plot_pdf(a_dur, color='k', linewidth=line_width_fit)

            a_dur1_pdf = pl.pdf(a_dur1, 10)
            a_dur2_pdf = pl.pdf(a_dur2, 10)
            BinCenters1 = (a_dur1_pdf[0][:-1]+a_dur1_pdf[0][1:])/2.
            BinCenters2 = (a_dur2_pdf[0][:-1]+a_dur2_pdf[0][1:])/2.
            x_max = a_dur1_pdf[0].max()
            interp1 = np.interp(np.arange(x_max), BinCenters1, a_dur1_pdf[1])
            interp2 = np.interp(np.arange(x_max), BinCenters2, a_dur2_pdf[1])
            fill_between(np.linspace(0,x_max,len(interp2)), interp2, interp1, facecolor='k', alpha=0.2)

        if regime == 'frozen':

            pl.plot_pdf(a_dur, color='cyan', linewidth=line_width)
            a_dur1_pdf = pl.pdf(a_dur1, 10)
            a_dur2_pdf = pl.pdf(a_dur2, 10)
            BinCenters1 = (a_dur1_pdf[0][:-1]+a_dur1_pdf[0][1:])/2.
            BinCenters2 = (a_dur2_pdf[0][:-1]+a_dur2_pdf[0][1:])/2.
            x_max = a_dur1_pdf[0].max()
            interp1 = np.interp(np.arange(x_max), BinCenters1, a_dur1_pdf[1])
            interp2 = np.interp(np.arange(x_max), BinCenters2, a_dur2_pdf[1])
            fill_between(np.linspace(0,x_max,len(interp2)), interp2, interp1, facecolor='cyan', alpha=0.2)

        subplot(222)
        if regime == 'non-frozen':
            pl.plot_pdf(a_area, color='k', linewidth=line_width_fit, \
                                                          label ='SORN')
            a_area1_pdf = pl.pdf(a_area1, 10)
            a_area2_pdf = pl.pdf(a_area2, 10)
            BinCenters1 = (a_area1_pdf[0][:-1]+a_area1_pdf[0][1:])/2.
            BinCenters2 = (a_area2_pdf[0][:-1]+a_area2_pdf[0][1:])/2.
            x_max = a_area1_pdf[0].max()
            interp1 = np.interp(np.arange(x_max), BinCenters1, a_area1_pdf[1])
            interp2 = np.interp(np.arange(x_max), BinCenters2, a_area2_pdf[1])
            fill_between(np.linspace(0,x_max,len(interp2)), interp2, interp1, facecolor='k', alpha=0.2)

        if regime == 'frozen':
            pl.plot_pdf(a_area, color='cyan', linewidth=line_width, \
                                                   label='Frozen (All)')
            a_area1_pdf = pl.pdf(a_area1, 10)
            a_area2_pdf = pl.pdf(a_area2, 10)
            BinCenters1 = (a_area1_pdf[0][:-1]+a_area1_pdf[0][1:])/2.
            BinCenters2 = (a_area2_pdf[0][:-1]+a_area2_pdf[0][1:])/2.
            x_max = a_area1_pdf[0].max()
            interp1 = np.interp(np.arange(x_max), BinCenters1, a_area1_pdf[1])
            interp2 = np.interp(np.arange(x_max), BinCenters2, a_area2_pdf[1])
            fill_between(np.linspace(0,x_max,len(interp2))[100:-140990], interp2[100:-140990], interp1[100:-140990], facecolor='cyan', alpha=0.2)



transient_steps = 0 # special way to plot an only initialized SORN
print 'Freezing plasticity in the beginning...'
exper_path =''
for regime in ['frozen']:

    print '\n', regime
    data_all = np.zeros((number_of_files, frozen_steps))

    for result_file in range(number_of_files):
        exper = 'result.h5'
        exper_path =  ''
        h5 = tables.openFile(os.path.join(exper_path,exper),'r')
        data = h5.root

        data_all[result_file] = np.around(data.activity[0] \
                             [frozen_steps+transient_steps:]*data.c.N_e)

        h5.close()

    a_dur, a_area = analysis.avalanches(data_all, 'N', '200', \
                                                        Threshold=THETA)


    subplot(221)
    pl.plot_pdf(a_dur, color='r', linewidth=line_width)


    subplot(222)
    pl.plot_pdf(a_area, color='r', linewidth=line_width, \
                                                 label='Random network')


# axis stuff
subplot(221)
xscale('log'); yscale('log')
xlim([1, 3000])
ylim([0.000001, 1])

# axis stuff
xlabel(r'$T$', fontsize=letter_size)
ylabel(r'$f(T)$', fontsize=letter_size)
fig_4a.xaxis.set_ticks_position('bottom')
fig_4a.yaxis.set_ticks_position('left')
tick_params(labelsize=letter_size)

fig_4a.spines['right'].set_visible(False)
fig_4a.spines['top'].set_visible(False)

# ticks name
xticks([1, 10, 100, 1000], \
       ['$10^{0}$', '$10^{1}$', '$10^{2}$', '$10^{3}$'])
yticks([1, 0.001, 0.000001], \
       ['$10^{0}$', '$10^{-3}$', '$10^{-6}$'])


subplot(222)
xscale('log'); yscale('log')
xlim([1, 20000])
ylim([0.0000001, 1])
xlabel(r'$S$', fontsize=letter_size)
ylabel(r'$f(S)$', fontsize=letter_size)
fig_4b.xaxis.set_ticks_position('bottom')
fig_4b.yaxis.set_ticks_position('left')
tick_params(labelsize=letter_size)
# ticks name
xticks([1, 100, 10000], \
       ['$10^{0}$', '$10^{2}$', '$10^{4}$'])
yticks([1, 0.001, 0.000001], \
       ['$10^{0}$', '$10^{-3}$', '$10^{-6}$'])

fig_4b.spines['right'].set_visible(False)
fig_4b.spines['top'].set_visible(False)

legend(loc=(0.4, 0.8), prop={'size':letter_size}, frameon=False)
########################################################################

########################################################################
# Fig. 4C, D

transient_steps = 2e6
frozen_steps = 2e6

for experiment_name in ['IPFrozen', 'AllButIPFrozen']:

    print experiment_name
    exper_path =''

    for regime in ['frozen']:

        print '\n', regime
        data_all = np.zeros((number_of_files, frozen_steps))

        for result_file in range(number_of_files):
            exper = 'result.h5'
            exper_path =  ''
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

        # shaded areas
        a_dur1, a_area1 = analysis.avalanches(data_all, 'N', '200',\
                                                      Theta_percent = 5)
        a_dur2, a_area2 = analysis.avalanches(data_all, 'N', '200',\
                                                     Theta_percent = 25)


        if experiment_name == 'IPFrozen':

            print 'Fig. 4C, D...\n'

            subplot(223)

            pl.plot_pdf(a_dur, color='r', linewidth=line_width)

            a_dur1_pdf = pl.pdf(a_dur1, 10)
            a_dur2_pdf = pl.pdf(a_dur2, 10)
            BinCenters1 = (a_dur1_pdf[0][:-1]+a_dur1_pdf[0][1:])/2.
            BinCenters2 = (a_dur2_pdf[0][:-1]+a_dur2_pdf[0][1:])/2.
            x_max = a_dur1_pdf[0].max()
            interp1 = np.interp(np.arange(x_max), BinCenters1, a_dur1_pdf[1])
            interp2 = np.interp(np.arange(x_max), BinCenters2, a_dur2_pdf[1])
            fill_between(np.linspace(0,x_max,len(interp2)), interp2, interp1, facecolor='r', alpha=0.2)



            subplot(224)
            label_exp = 'Frozen (IP)'
            pl.plot_pdf(a_area, linewidth=line_width, color='r', \
                                                        label=label_exp)

            a_area1_pdf = pl.pdf(a_area1, 10)
            a_area2_pdf = pl.pdf(a_area2, 10)
            BinCenters1 = (a_area1_pdf[0][:-1]+a_area1_pdf[0][1:])/2.
            BinCenters2 = (a_area2_pdf[0][:-1]+a_area2_pdf[0][1:])/2.
            x_max = a_area1_pdf[0].max()
            interp1 = np.interp(np.arange(x_max), BinCenters1, a_area1_pdf[1])
            interp2 = np.interp(np.arange(x_max), BinCenters2, a_area2_pdf[1])
            fill_between(np.linspace(0,x_max,len(interp2)), interp2, interp1, facecolor='r', alpha=0.2)


        if experiment_name == 'AllButIPFrozen':

            subplot(223)
            pl.plot_pdf(a_dur, color='cyan', linewidth=line_width)

            a_dur1_pdf = pl.pdf(a_dur1, 10)
            a_dur2_pdf = pl.pdf(a_dur2, 10)
            BinCenters1 = (a_dur1_pdf[0][:-1]+a_dur1_pdf[0][1:])/2.
            BinCenters2 = (a_dur2_pdf[0][:-1]+a_dur2_pdf[0][1:])/2.
            x_max = a_dur1_pdf[0].max()
            interp1 = np.interp(np.arange(x_max), BinCenters1, a_dur1_pdf[1])
            interp2 = np.interp(np.arange(x_max), BinCenters2, a_dur2_pdf[1])
            fill_between(np.linspace(0,x_max,len(interp2)), interp2, interp1, facecolor='cyan', alpha=0.2)

            subplot(224)
            label_exp = 'Frozen (All but IP)'
            pl.plot_pdf(a_area, color='cyan', linewidth=line_width, \
                                                        label=label_exp)

            a_area1_pdf = pl.pdf(a_area1, 10)
            a_area2_pdf = pl.pdf(a_area2, 10)
            BinCenters1 = (a_area1_pdf[0][:-1]+a_area1_pdf[0][1:])/2.
            BinCenters2 = (a_area2_pdf[0][:-1]+a_area2_pdf[0][1:])/2.
            x_max = a_area1_pdf[0].max()
            interp1 = np.interp(np.arange(x_max), BinCenters1, a_area1_pdf[1])
            interp2 = np.interp(np.arange(x_max), BinCenters2, a_area2_pdf[1])
            fill_between(np.linspace(0,x_max,len(interp2)), interp2, interp1, facecolor='cyan', alpha=0.2)



# axis stuff
subplot(223)
xscale('log'); yscale('log')
xlim([1, 3000])
ylim([0.000001, 1])
xlabel(r'$T$', fontsize=letter_size)
ylabel(r'$f(T)$', fontsize=letter_size)
fig_4c.xaxis.set_ticks_position('bottom')
fig_4c.yaxis.set_ticks_position('left')
tick_params(labelsize=letter_size)
# ticks name
xticks([1, 10, 100, 1000], \
       ['$10^{0}$', '$10^{1}$', '$10^{2}$', '$10^{3}$'])
yticks([1, 0.001, 0.000001], \
       ['$10^{0}$', '$10^{-3}$', '$10^{-6}$'])



subplot(224)
xscale('log'); yscale('log')
xlim([1, 20000])
ylim([0.0000001, 1])
xlabel(r'$S$', fontsize=letter_size)
ylabel(r'$f(S)$', fontsize=letter_size)
# ticks name
xticks([1, 100, 10000], \
       ['$10^{0}$', '$10^{2}$', '$10^{4}$'])
yticks([1, 0.001, 0.000001], \
       ['$10^{0}$', '$10^{-3}$', '$10^{-6}$'])

fig_4c.spines['right'].set_visible(False)
fig_4c.spines['top'].set_visible(False)
fig_4d.spines['right'].set_visible(False)
fig_4d.spines['top'].set_visible(False)

fig_4d.xaxis.set_ticks_position('bottom')
fig_4d.yaxis.set_ticks_position('left')
tick_params(labelsize=letter_size)

legend(loc=(0.3, 0.8), prop={'size':letter_size}, frameon=False)
########################################################################
fig_4a.annotate('A', xy=subplot_letter, xycoords='axes fraction', \
                fontsize=letter_size_panel ,  fontweight='bold', \
                horizontalalignment='right', verticalalignment='bottom')
fig_4b.annotate('B', xy=subplot_letter, xycoords='axes fraction', \
                fontsize=letter_size_panel ,  fontweight='bold', \
                horizontalalignment='right', verticalalignment='bottom')
fig_4c.annotate('C', xy=subplot_letter, xycoords='axes fraction', \
                fontsize=letter_size_panel ,  fontweight='bold', \
                horizontalalignment='right', verticalalignment='bottom')
fig_4d.annotate('D', xy=subplot_letter, xycoords='axes fraction', \
                fontsize=letter_size_panel ,  fontweight='bold', \
                horizontalalignment='right', verticalalignment='bottom')
fig.subplots_adjust(wspace=.3)
fig.subplots_adjust(hspace=.3)

print 'Saving figures...',
result_path = '../../plots/'
result_name_png = 'Fig4.pdf'
savefig(os.path.join(result_path, result_name_png), format = 'pdf')
print 'done\n\n'

show()
