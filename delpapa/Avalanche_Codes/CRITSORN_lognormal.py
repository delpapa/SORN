####
# lognormal distributions plots 
####

from pylab import * 
from scipy.optimize import curve_fit # to fit the log normal curves

import tables
import os
from tempfile import TemporaryFile

import data_analysis as analysis                                   

figure(1, figsize=(12,12))

number_of_files = 10

# figure parameters
letter_size = 25
line_width = 4.0
points_size = 13
subplot_letter = (0.95, 0.9)

########################################################################
experiment_name = 'Lognormal_Stuff'
print experiment_name

for sigma in xrange(number_of_files):

    print sigma
        
    exper = 'result.h5'
    exper_path =  '../Avalanche_Experiments/Lognormal_Stuff/Average/'\
         + sigma + '/common/'
    h5 = tables.openFile(os.path.join(exper_path,exper),'r')
    data = h5.root
    logweight = data.endweight[0][data.endweight[0]>0]   
    import ipdb; ipdb.set_trace()
    h5.close()
    
    for bining_size in [50]: # [5, 10, 20, 50, 100]:    
        logbins = logspace(-2,0, bining_size)
        (y,_) = histogram(logweight, bins=logbins)
    
        x = logbins[:-1]+(logbins[0]+logbins[1])/2.0    #fit data to lognormal
        semilogx(x,y,'.')

    
        # Do the fitting
        def lognormal(x,mue,var,scale): 
            return scale * (exp(- ((log(x)-mue)*(log(x)-mue)) /\
                                          (2*var)) / (x*sqrt(2*pi*var)))
                                          
        y=y.astype('float')/y.sum()  # normalization
                                    
        #popt, pcov = curve_fit(lognormal, x, y)
        #curve_x = logspace(-2,0,100)
        #fitted_y = lognormal(curve_x,*popt)
    

        # if bining_size == 10:
            # plot(curve_x,fitted_y, 'k', linewidth = line_width, label = 'lognormal for 10')
        plot(x,y, '-o', markersize=points_size, linewidth = line_width, label=sigma)

            
xlim([0.02, 1])
ylim([0, 0.3])

xlabel('Synaptic Weight', fontsize=letter_size)
ylabel('Frequency', fontsize=letter_size)
legend(loc='best', prop={'size':letter_size}, frameon=False)
#~ setp(ax, xticks=[0.01, 0.1, 1], xticklabels=['$10^{-2}$', '$10^{-2}$', '$10^0$'], yticks=[4500, 5000], yticklabels=['$4500$', '$5000$'])
#~ setp(ax2, xticks=[0.01, 0.1, 1], xticklabels=['$10^{-2}$', '$10^{-2}$', '$10^0$'], yticks=[0, 500, 1000], yticklabels=['$0$', '$500$', '$1000$'])

tick_params(axis='both', which='major', labelsize=letter_size)

########################################################################

tight_layout()

print 'Saving figure...',		
result_path = '../Avalanche_Results/'
result_name_png = 'CRITSORN_Lognormal_etaSTDP.png'

savefig(os.path.join(result_path, result_name_png), format = 'png')
print 'done'
show()
