####
# Script for the sixth paper figure
# Includes: noise-ABCD-noise
####

from pylab import * 

import tables
import os
from tempfile import TemporaryFile
from matplotlib.colors import LogNorm

import data_analysis as analysis


number_of_files = 1
########################################################################
# Counting Task - Performance                                          #
########################################################################
print '\nCalculating performance for the Counting Task...'
            
final_performance_mean = []
final_performance_std = []
final_performance_5p = []
final_performance_95p = []
final_sequence_lengh = [6]
for experiment_folder in ['n6_andrea']:

    print experiment_folder

    partial_performance = np.zeros(number_of_files)
    
    for file_number in range(number_of_files):
        
        # read data files 
        exper = 'result.h5'
        exper_path =  '../../../Avalanche_Experiments/Counting_task/' + \
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
    final_performance_5p.append(np.percentile(partial_performance, 16))
    final_performance_95p.append(np.percentile(partial_performance, 84))

figure()
imshow(W_trained, interpolation='none', cmap='gray_r', aspect='auto', norm = LogNorm())  
# colorbar(norm = LogNorm())       
show()
