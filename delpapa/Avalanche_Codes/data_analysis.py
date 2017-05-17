########################################################################
# This script contains all the data analysis functions                 #
#                                                                      #
# Bruno Del Papa, 23.07.2015                                           #
########################################################################

from __future__ import division
from pylab import *
import scipy, scipy.stats

import tables
import os
from tempfile import TemporaryFile

### log binning function
# calculate the log binning for a single curve
# return the (x,y) coordinates of the binned points and the (x,y) data 
# for the non-binned (but normalized) data
# - data MUST be an array of integers!!
def log_binning(data, variable, value,
		        bin_size = 0.1): # bin size for the binned data

					
    # Way to plot the frequency x avalanche size. 
    # stackoverflow.com/questions/10119441/numpy-bincount-with-floats
    x_points, inverse = unique(data, return_inverse=True)
    y_freq = bincount(inverse)
    y_freqs_norm = y_freq / float(y_freq.sum()) # normalization
    
    bin_x = []; bin_y = []
    for i in arange(log10(x_points.min()), log10(x_points.max()), \
                                                              bin_size):
        elements = where(np.logical_and(i <= log10(x_points), \
                                          log10(x_points) < i+bin_size))
        if x_points[elements].size != 0:
            bin_x.append(x_points[elements].mean())
            bin_y.append(y_freqs_norm[elements].mean())
    
    # Here bin_x, bin_y are the ORDERED coordinates of the binned data        
    bin_x = asarray(bin_x)
    bin_y = asarray(bin_y)

    return bin_x, bin_y, x_points, y_freqs_norm

### peak's avalanche (fixed Theta)
# return 2 arrays with duration and area of all avalanches, in the 
# same order (although non-ordered by anything except order of 
# apearence)
# - activity is a matrix for all the experiments (one in each row)
# - Theta, if not given, is 15%, about half of the average activity
# (Poli et al.). Otherwise is a percentile
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

### plot avalanches from files
# This function do not save the Figures, which must be saved in the 
# main code
def plot_aval_files(variable, v):
	
    # plot parameters
    dur_max = 1e4
    area_max = 1e5
	
    if variable == 'p':
        result_path = '../Avalanche_Results/RandomSpikes/'+ variable + v
    else:                            
        result_path = '../Avalanche_Results/Aval_Comp/' + variable + v
        	
    bin_d_x = fromfile(result_path +'_dur_x.txt', sep = ' ')
    bin_d_y = fromfile(result_path +'_dur_y.txt', sep = ' ')
    bin_a_x = fromfile(result_path +'_area_x.txt', sep = ' ')
    bin_a_y = fromfile(result_path +'_area_y.txt', sep = ' ')	
	
    ### duration Figure
    figure(1, figsize=(12,5.5))	
    subplot(121)
    xlim(1, dur_max)
    dur_range = arange(dur_max)[1:]
    
    if variable == 'N':
        plot(bin_d_x, bin_d_y, linewidth = 2.0, \
                                            label=variable+r'$^E$ = '+v)
                                            
    elif variable == 'Hip':                                      
        plot(bin_d_x, bin_d_y, linewidth = 2.0, \
                                                  label=r'$H_{IP}$ ='+v)                                        
                                            
    elif variable == 'eta':
        plot(bin_d_x, bin_d_y, linewidth = 2.0, \
                                         label=r'$ \eta _{IP}$'+' = '+v)
                                            
    elif variable == 'sigma':
        plot(bin_d_x, bin_d_y, linewidth = 2.0, \
                                            label=variable+r'$^E$ = '+v)          
        
    xscale('log'); yscale('log')
    xlabel('T', fontsize=15);
    ylabel('f(T)', fontsize=15) 
    tick_params(axis='both', which='major', labelsize=15)
    
    ### area Figure
    subplot(122)
    xlim(1, area_max)
    area_range = arange(area_max)[1:]	
    if variable == 'N':                                      
        plot(bin_a_x, bin_a_y, linewidth = 2.0, \
                                            label=variable+r'$^E$ = '+v)
                                            
    elif variable == 'Hip':                                      
        plot(bin_a_x, bin_a_y, linewidth = 2.0, \
                                                  label=r'$H_{IP}$ ='+v)
                                            
    elif variable == 'eta':
        plot(bin_a_x, bin_a_y, linewidth = 2.0, \
                                            label=r'$ \eta _{IP}$ = '+v)
                                         
    elif variable == 'sigma':
        plot(bin_a_x, bin_a_y, 'k', linewidth = 2.0, \
                                           label=r'$\sigma^2_\xi$ = '+v)

    xscale('log'); yscale('log')	 
    xlabel('S', fontsize = 15)
    ylabel('f(S)', fontsize=15)
    legend(loc='best', prop={'size':15})  
    tick_params(axis='both', which='major', labelsize=15)			

### plot activity from files
# This function do not save the Figures, which must be saved in the 
# main code. 
def plot_act_files(variable, v):
    
    if variable == 'p':
        result_path = '../Avalanche_Results/RandomSpikes/'+ variable + v
    else:                            
        result_path = '../Avalanche_Results/Aval_Comp/' + variable + v
        	
    dist_mean = fromfile(result_path +'_activity.txt', sep = ' ')
    dist_std = fromfile(result_path +'_activity_std.txt', sep = ' ')
    
    # activity distribution figure
    figure(2)
    
    if variable == 'N':
        plot(dist_mean, linewidth = 2.0, label=variable+r'$^E$ = '+v)
    elif variable == 'Hip':
        plot(dist_mean, linewidth = 2.0, label=r'$H_{IP}$ = '+v)
    elif variable == 'eta':
        plot(dist_mean, linewidth = 2.0, label=r'$ \eta _{IP}$ = '+v)
    elif variable == 'sigma':
        plot(dist_mean, linewidth = 2.0, label=r'$\sigma^2_\xi$ = '+v)

    xlabel('activity', fontsize=15); 
    ylabel('activity distribution', fontsize=15); 
    legend(loc='best', prop={'size':15})
    tick_params(axis='both', which='major', labelsize=15)
    xlim([0,80])
        
### plot the entropy for ONE set (variable, value) - 03.12.14
### activity is the activity data for all files
### save 3 txt files with the xrange and the dur and area entropy
def entropy_allTheta(activity, variable, value):
	
	entropy_area = []; entropy_dur = []
	ent_max_range = int(float(value))+1
	x_range = np.arange(ent_max_range+1)
	for theta in x_range:
		sys.stdout.write('\rcalculating avalanches %d%%' \
											%(100*theta/ent_max_range)),
		sys.stdout.flush()
		a_duration_all = []; a_area_all = []
		for data_file in range(len(activity)):
			avalan_stable = activity[data_file] - theta
			size, area = 0, 0
			for i in range (len(avalan_stable)):
				if avalan_stable[i] > 0:
					size += 1
					area += int(avalan_stable[i]) # to avoid errors
				elif size != 0:
					a_duration_all.append(size)
					a_area_all.append(area)
					size, area = 0, 0
		a_duration = np.asarray(a_duration_all)
		a_area = np.asarray(a_area_all)	
		
		# calculate the normalized sorted avalanche distribution
		# TODO there must be a better way 
		dur_points, inverse = np.unique(a_duration, return_inverse=True)
		dur_freq = np.bincount(inverse) 
		dur_freq_norm = dur_freq / float(dur_freq.sum())			
		area_points, inverse = np.unique(a_area, return_inverse=True)
		area_freq = np.bincount(inverse)
		area_freq_norm = area_freq / float(area_freq.sum())	
		entropy_dur.append( \
						   (-dur_freq_norm*np.log(dur_freq_norm)).sum())
		entropy_area.append( \
						 (-area_freq_norm*np.log(area_freq_norm)).sum())
	entropy_dur = np.asarray(entropy_dur)
	entropy_area = np.asarray(entropy_area)
	
	result_path = '../Avalanche_Results/Entropy/' + variable + value
	x_range.tofile(result_path +'_xrange.txt', sep = ' ')
	entropy_dur.tofile(result_path +'_entr_dur.txt', sep = ' ')
	entropy_area.tofile(result_path +'_entr_area.txt', sep = ' ')
	print '...done'
	return x_range, entropy_dur, entropy_area

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

### plot the best power-law fit (MLE - [Clauset et al. 2009])
# Use as input a set of a discreet distribution (non-ordered)
# Return exponent alpha and alpha error
def power_law_disc(data, x_min, plot_fit = True):
    
    data_effec = data[x_min <= data]
    n = len(data_effec)
    alpha = 1 + n / (sum(log(data_effec/(x_min-0.5))))
    sigma_alpha = (alpha-1)/sqrt(n)
    
    x_range = arange(1, data_effec.max())
    if plot_fit:
        plot(x_range, x_range**-alpha, '-b', label = r'$ \alpha$ =' + \
                                        str(around(-alpha, decimals=3)))
    
    return alpha, sigma_alpha
    
def fractal_plot(activity, Theta = False):
	
	if Theta == False:
		Theta = int(activity.mean()/2)
	avalan_stable = activity - Theta
	a_duration = []; a_area = []; a_curves = []
	size, area, curve = 0, 0, []
	for i in range (len(avalan_stable)):
		if avalan_stable[i] > 0:
			size += 1
			area += int(avalan_stable[i]) # to avoid rouding errors
			curve.append(int(avalan_stable[i]))
		elif size != 0:
			a_duration.append(size)
			a_area.append(area)
			curve.insert(0,0); curve.append(0)
			a_curves.append(curve)
			size, area, curve = 0, 0, []
	a_curves.sort(key = len)
	import ipdb; ipdb.set_trace()
	
		

	figure()
	for c in a_curves:
		c = np.asarray(c)
		c_l = float(len(c))
		c_max = c.max()
		x_bin = []; y_bin = []
		for i, j in enumerate(c):
			x_bin.append(i/(c_l-1)) 
			y_bin.append(j/c_max)
		if i%50 == 0:	
			plot(x_bin, y_bin)
	show()

