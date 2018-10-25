"""Fits Adaptive Integrate-and-Fire (AdEx) model to all the data files in the
current map and saves the results in data files for inh/exc adex
- Matlab file 'export_data.m" is used to convert the matlab struct to
three datafiles containing spikeindices, input current and membrane
potential and to give the files the correct format.
- Only the data containing more than min_spikes size are used
and 2/3 of the data (or 50 spikes if 2/3rd is >50) is used for
optimization with brians modelfitting function.
- numpy 1.10.4 is used with python 2.7, a newer version of nummpy
will result in errors when using Brians modelfitting function.
- Library of functions in 'optimize_vpeak.py'

TODO:
Set output_file_adex at new filename

SECOND VERSION, Vpeak added
26-05-2018, Linda Wouters
"""


# Importing optimization functions from library 'optimize_vpeak',
# since modelfitting gives the exact same incorrect value for
# best_fit if called multiple times within the same script
import optimize_vpeak

import numpy as np
import sys
import glob

# Hardcoded, dependend on type of data
SAMPLE_RATE = 20.


# Still gives good fits with short run time
popsize = 5000
# No difference in gamma between 10, 20, 50 or 150 maxiters
# No difference between 20 and 50 for izhik
maxiter = 20

N = 1
delta = 4
cpus = 6
min_spikes = 40
# TODO: Check output file names
output_file_adex = "optim_multi_adex_2605_0.csv"
output_file_izhik = "optim_multi_izhik_2605_0.csv"

# if_name_==__main__ to avoid crashing 
# http://www.briansimulator.org/docs/modelfitting.html
if __name__=='__main__':

    for i in range(N):
        # Iterate over all the files in the current directory ending with 'analyzed_spikeindices.csv"
        # Code from: https://stackoverflow.com/questions/18262293/how-to-open-every-file-in-a-folder
        # answered Aug 15 '13 at 21:38 by Viktor Kerkez
        for spike_file in glob.glob('*analyzed_spikeindices.csv'):
            # Use try/except to carry on to the next file if any error occurs
            try:
                print spike_file
                # Get current_file
                current_file = spike_file.replace("spikeindices", "input")

                # Get spike times and current input
                spikes, current, __ = optimize_vpeak.read_data(spike_file,
                                                current_file,
                                                SAMPLE_RATE)

                # Downsample to half the sample rate
                # Code from: https://stackoverflow.com/questions/34231244/downsampling-a-2d-numpy-array-in-python
                # answered Dec 11 '15 at 21:02 by Bart 
                # from 7200000 to 3600000 length
                current = current[::2]

                # Set downsampled sample_rate
                # assert isinstance(SAMPLE_RATE, float)
                data_dt = (2/float(SAMPLE_RATE))

                # Extract neuron_type from file name
                neuron_type = spike_file[3:6]
                # print "neuron_type:", neuron_type
                if (neuron_type != "exc" and neuron_type != "inh"):
                    raise ValueError("Neurontype is not 'exc' or 'inh'.")

                # Optimize if there are more than min_spikes spikes
                if len(spikes) > min_spikes and (neuron_type == "exc" or neuron_type =="inh"):
                    if len(spikes)*2/3 > 51:
                        # Set fitrange to 50th spike
                        fit_range = spikes[50] / data_dt
                        print("50", fit_range)
                    else:
                        # Set fitrange to the 2/3rd spikelength-spike
                        fit_range = spikes[int(len(spikes)*2/3)] / data_dt
                        print("2/3", fit_range)
                
                    # Shorten data to speed up optimization
                    spikes = np.array([x for x in spikes if x < fit_range*data_dt])
                    print len(spikes), spikes
                    current = current[:int(fit_range)]

                    # Get results
                    results = optimize_vpeak.modelfit("adex", spikes, current, neuron_type,
                                                data_dt, popsize, maxiter, delta, cpus,
                                                spike_file, current_file, output_file_adex)
                    print results

                    # results = optimize_vpeak.modelfit("izhik", spikes, current, neuron_type,
                    #                             data_dt, popsize, maxiter, delta, cpus,
                    #                             spike_file, current_file, output_file_izhik)
                    # print results

                # If spikes length was < min_spikes, write None for all values for this file
                else:
                    optimize_vpeak.write_results("adex", None, output_file_adex, spike_file)
                    # optimize_vpeak.write_results("izhik", None, output_file_izhik, spike_file)

            # Except code from: https://wiki.python.org/moin/HandlingExceptions#CA-1dd2cf7976df10fc545caa391478799340411158_1
            # paragraph "General Error Catching"
            except KeyboardInterrupt as error:
                sys.exit()
            except:
                print "Error: %s" % sys.exc_info()[1]
                # print sys.exc_info()
                sys.exc_clear()



