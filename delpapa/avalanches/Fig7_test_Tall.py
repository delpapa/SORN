####
# Script for the sixth paper figure
# Includes: noise-ABCD-noise
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

### figure parameters
width  =  7
height = 8
fig_7 = figure(1, figsize=(width, height))

letter_size = 16
letter_size_panel = 16
line_width = 1.5
line_width_fit = 2.0
subplot_letter = (-0.1, 0.9)
subplot_letter1 = (-0.3 , 0.95)

number_of_files = 10

########################################################################
# Counting Task - Performance                                          #
########################################################################
print '\nCalculating performance for the Counting Task...'

title('Counting Task')           
            
final_performance_mean = []
final_performance_std = []
final_performance_5p = []
final_performance_95p = []
final_sequence_lengh = [4, 6, 8, 14, 20]
for experiment_folder in ['n4_andrea', 'n6_andrea', 'n8_andrea', \
                           'n14_andrea', 'n20_andrea']:

    print experiment_folder

    partial_performance = np.zeros(number_of_files)
    
    for file_number in range(number_of_files):
        
        # read data files 
        exper = 'result.h5'
        exper_path =  '../Avalanche_Experiments/Counting_task/' + \
                 experiment_folder + '/' + str(file_number+1) + '/common/'
        h5 = tables.openFile(os.path.join(exper_path,exper),'r')
        data = h5.root

        # create training and test arrays 
        train_steps = data.c.steps_readouttrain[0]
        test_steps = data.c.steps_readouttest[0]

        # the letter to be presented at the time step
        y_train = data.countingletter[0][:train_steps]  
        y_test = data.countingletter[0][train_steps:train_steps+test_steps]

        # divide readout into n different units
        y_read_train = np.zeros((6, len(y_train)))
        y_read_test = np.zeros((6, len(y_test)))
        for i, y in enumerate(y_train):
            y_read_train[y, i] = 1
        for i, y in enumerate(y_test):
            y_read_test[y, i] = 1
        target = np.argmax(y_read_test, axis=0) # target for the training data

        # internal state before letter presentatio
        X_train = (data.countingactivity[0][:,:train_steps] >= 0) + 0. 
        X_test = (data.countingactivity[0][:,train_steps:train_steps+test_steps] >= 0) + 0.  

        h5.close()

        # Readout training
        X_train_pinv = np.linalg.pinv(X_train) # MP pseudo-inverse
        W_trained = np.dot(y_read_train, X_train_pinv) # least squares weights

        # Network prediction with trained weights
        y_predicted = np.dot(W_trained, X_test)

        # performance
        prediction = np.argmax(y_predicted, axis=0)          
        perf_all = (prediction == target).sum()/float(len(y_test))

        # reduced performance (removing the prediction of the first letter)
        except_first = np.where(np.logical_or(\
                            np.logical_or(y_test == 1, y_test == 2),\
                            np.logical_or(y_test == 4, y_test == 5)))[0]
        y_test_red = y_test[except_first]
        y_pred_red = prediction[except_first]
        perf_red = (y_test_red==y_pred_red).sum()/float(len(y_pred_red))
        
        partial_performance[file_number] = perf_red
        
    final_performance_mean.append(partial_performance.mean())
    final_performance_std.append(partial_performance.std()) 
    final_performance_5p.append(np.percentile(partial_performance, 25))
    final_performance_95p.append(np.percentile(partial_performance, 75))

subplot(412)     
plot(final_sequence_lengh, final_performance_mean, '-k', \
                                               label=r'Original SORN')
low_err = np.array(final_performance_mean)-np.array(final_performance_5p)
hig_err = np.array(final_performance_95p)-np.array(final_performance_mean)                                               
errorbar(final_sequence_lengh, final_performance_mean, \
                                  yerr=[low_err, hig_err], color='k')
                                  
                                  
final_performance_mean = []
final_performance_std = []
final_performance_5p = []
final_performance_95p = []
final_sequence_lengh = [4, 6, 8, 14, 20]
for experiment_folder in ['n4_Tall', 'n6_Tall', 'n8_Tall', \
                           'n14_Tall', 'n20_Tall']:

    print experiment_folder

    partial_performance = np.zeros(number_of_files)
    
    for file_number in range(number_of_files):
        
        # read data files 
        exper = 'result.h5'
        exper_path =  '../Avalanche_Experiments/Counting_task/' + \
                 experiment_folder + '/' + str(file_number+1) + '/common/'
        h5 = tables.openFile(os.path.join(exper_path,exper),'r')
        data = h5.root

        # create training and test arrays 
        train_steps = data.c.steps_readouttrain[0]
        test_steps = data.c.steps_readouttest[0]

        # the letter to be presented at the time step
        y_train = data.countingletter[0][:train_steps]  
        y_test = data.countingletter[0][train_steps:train_steps+test_steps]

        # divide readout into n different units
        y_read_train = np.zeros((6, len(y_train)))
        y_read_test = np.zeros((6, len(y_test)))
        for i, y in enumerate(y_train):
            y_read_train[y, i] = 1
        for i, y in enumerate(y_test):
            y_read_test[y, i] = 1
        target = np.argmax(y_read_test, axis=0) # target for the training data

        # internal state before letter presentatio
        X_train = (data.countingactivity[0][:,:train_steps] >= 0) + 0. 
        X_test = (data.countingactivity[0][:,train_steps:train_steps+test_steps] >= 0) + 0.  

        h5.close()

        # Readout training
        X_train_pinv = np.linalg.pinv(X_train) # MP pseudo-inverse
        W_trained = np.dot(y_read_train, X_train_pinv) # least squares weights

        # Network prediction with trained weights
        y_predicted = np.dot(W_trained, X_test)

        # performance
        prediction = np.argmax(y_predicted, axis=0)          
        perf_all = (prediction == target).sum()/float(len(y_test))

        # reduced performance (removing the prediction of the first letter)
        except_first = np.where(np.logical_or(\
                            np.logical_or(y_test == 1, y_test == 2),\
                            np.logical_or(y_test == 4, y_test == 5)))[0]
        y_test_red = y_test[except_first]
        y_pred_red = prediction[except_first]
        perf_red = (y_test_red==y_pred_red).sum()/float(len(y_pred_red))
        
        partial_performance[file_number] = perf_red
        
    final_performance_mean.append(partial_performance.mean())
    final_performance_std.append(partial_performance.std()) 
    final_performance_5p.append(np.percentile(partial_performance, 25))
    final_performance_95p.append(np.percentile(partial_performance, 75))

subplot(412)     
plot(final_sequence_lengh, final_performance_mean, '-r', \
                                               label=r'Same T update')
low_err = np.array(final_performance_mean)-np.array(final_performance_5p)
hig_err = np.array(final_performance_95p)-np.array(final_performance_mean)                                               
errorbar(final_sequence_lengh, final_performance_mean, \
                                  yerr=[low_err, hig_err], color='r')                                  
                                  

ylabel('Performance', fontsize=letter_size)
xlabel('n', fontsize=letter_size)
legend(loc='best', prop={'size':letter_size}, frameon=False)
xlim([0, 25])
ylim([0.4, 1.1])
tick_params(axis='both', which='major', labelsize=letter_size)
xticks([0, 5, 10, 15, 20, 25], \
     ['$0$', '$5$', '$10$', '$15$', '$20$', '$25$'])
yticks([0.4, 0.6, 0.8, 1.0],\
            ['$0.4$', '$0.6$', '$0.8$', '$1.0$'])  
            
figure(2)
exper = 'result.h5'
exper_path =  '../Avalanche_Experiments/Counting_task/' + \
                 'n20_andrea' + '/' + str(2) + '/common/'
h5 = tables.openFile(os.path.join(exper_path,exper),'r')
data = h5.root
raster = data.Spikes[0][:, -1000:]
for (i,sp) in enumerate(raster):
    s_train = where(sp == 1)[0]
    if s_train != []:
        vlines(s_train, i + 0.5, i + 1.5)
        hold('on')
ylabel('Low', fontsize=letter_size)
ylim([0, 200])
tick_params(axis='both', which='major', labelsize=letter_size) 
xticks([])
yticks([])                      
            
########################################################################
print 'Saving figures...',		
result_path = '../Avalanche_Results/'
result_name_png = 'Fig7_test_Tall.pdf'
savefig(os.path.join(result_path, result_name_png), format = 'pdf')

show()
