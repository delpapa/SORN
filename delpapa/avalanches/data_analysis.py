########################################################################
# This script contains all the data analysis functions                 #
########################################################################

from __future__ import division
from pylab import *
import scipy, scipy.stats

import tables
import os
from tempfile import TemporaryFile

def avalanches(activity, variable, value,
		       Threshold = 'percent', Theta_percent = 15,
               Transient = 0, fullS = False, binsize = False):

    # Threhold indicates which kind of activity treshold are we using
    # If FALSE, Theta_percent percentile is used
    # if 'half', 'half' the mean activity is used
    # if '1std', mean activity minus one std is used

    a_duration_all = []; a_area_all = []


    # Theta is by default given by 25 percentile
    # Theta MUST be a int
    if Threshold == 'percent':
        Theta = percentile(activity, Theta_percent)
    elif Threshold == 'half':
        Theta = int(activity.mean()/2.)
    elif Threshold == '1std':
        Theta = int(activity.mean()-activity.std())
    else:
        Theta = Threshold

    for data_file in range(len(activity)):

        sys.stdout.write('\rcalculating avalanches %d%%' \
                                    %(100*(data_file+1)/len(activity))),
        sys.stdout.flush()

        # make prettier
        if binsize is not False and binsize != 1:

            avalan_stable = activity[data_file]

            if binsize == 2:
                avalan_stable = avalan_stable[::2] + avalan_stable[1::2]
                avalan_stable /= 2.

            if binsize == 5:
                avalan_stable = avalan_stable[::5] + avalan_stable[1::5] +\
                                avalan_stable[2::5] + avalan_stable[3::5] +\
                                avalan_stable[4::5]
                avalan_stable /= 5.
            if binsize == 10:
                avalan_stable = avalan_stable[::10] + avalan_stable[1::10] +\
                                avalan_stable[2::10] + avalan_stable[3::10] +\
                                avalan_stable[4::10] + avalan_stable[5::10] +\
                                avalan_stable[6::10] + avalan_stable[7::10] +\
                                avalan_stable[8::10] + avalan_stable[9::10]
                avalan_stable /= 10.

            avalan_stable = np.floor(avalan_stable - Theta)

        else:
            # to avoid empty array error
            if len(activity.shape) > 1:
                avalan_stable = activity[data_file] - Theta
            else:
                avalan_stable = activity - Theta

        size, area = 0, 0
        ### TODO: make it prettier - this should include only the avalanches
        ### before the transient
        transient_end = 0
        new_range = len(avalan_stable)
        if Transient is not 0:
            for i in range(Transient, new_range):
                if avalan_stable[i] == 0:
                    transient_end = i + 1
                    break
            new_range = transient_end
        ####

        for i in range(new_range):
            if avalan_stable[i] > 0:
                size += 1
                if not fullS:
                    area += int(avalan_stable[i])
                else:
                    area += int(activity[data_file][i])
            elif size != 0:
                a_duration_all.append(size)
                a_area_all.append(area)
                size, area = 0, 0


    # convert to np.array cause it is easier for the other functions
    a_duration_all = asarray(a_duration_all)
    a_area_all = asarray(a_area_all)
    print '...done'
    return a_duration_all, a_area_all

### distribution of the total activity
# Return the average activity (as a array), and std
def mean_activity(activity, variable, value):

	distribution = zeros((len(activity), activity.max()+1))
	for data_file in range(len(activity)):

        # print the % in the terminal
		sys.stdout.write('\rcalculating activity %d%%' \
									%(100*(data_file+1)/len(activity))),
		sys.stdout.flush()

		total_steps = activity[data_file].size
		for i in range(total_steps):
			distribution[data_file, activity[data_file, i]] += 1
		distribution[data_file, :] /= distribution[data_file, :].sum()

	dist_mean = distribution.mean(0)
	dist_std = distribution.std(0)
	print '...done'

	return dist_mean, dist_std

### calculate the size average as a function of the duration
# receives the non-sorted arrays with measures of size and duration
# returns two non-sorted arrays containing the duration and average
# avalanche size.
def area_X_duration(a_dur, a_area):

    S_avg = []
    T_avg = []

    for i in range(len(a_dur)):
        duration = a_dur[i]
        if duration not in T_avg:
            T_avg.append(duration)
            S_avg.append(a_area[np.where(a_dur==duration)].mean())

    T_avg=asarray(T_avg)
    S_avg=asarray(S_avg)

    return T_avg, S_avg
