"""Runs the network with three neurons per population.
Uses 'network_library_threepop.py' as library for the functions

TODO:
-Change the output_file_ctrl and output_file_depr
-Set N_iter for the number of simulations

Usage: 
python network_threepop.py control
or
python network_threepop.py deprived

06-2018, Linda Wouters"""

from brian import *
import network_library_threepop as lib
import csv
import os
import peakutils

import matplotlib.pyplot as plt


output_file_ctrl = "network3_181023_ctrl_avg1_0_sp"
output_file_depr = "network3_181023_depr_avg1_0_sp"

N_iter = 5*10**4

# +-.05mV, else the first peak is not detected
peak_thres = .00005
sampling_rate = .1*ms/second
# number of indices around the spike as range
spike_thres = 2

# No averaging
avg_N = 1

max_gsyn = 7.5*nsiemens
plot = False

if __name__=='__main__':

    # Get equations and parameters for right model
    if len(sys.argv) != 2:
        print "expected command argument control/deprived"
        sys.exit()
    elif sys.argv[1] == "control":
        deprived = False
        output_file = output_file_ctrl
    elif sys.argv[1] == "deprived":
        deprived = True
        output_file = output_file_depr
    else:
        print "expected command argument control/deprived"
        sys.exit()
    # Add counter to filename, so multiple files can be created
    # in case of a memoryerror
    file_counter = 0
    output_file_name = output_file+ '_p' + str(file_counter) + '.csv'

    # Initialize equations, parameters and variables
    eqs, reset_eqs, thres_eqs, syn_model, syn_eqs, syn_onpre, params_e,\
        params_i, variables, Vinit, dt, N = lib.get_initial_vals()
    # Set headers for the datafile
    headers = ['gsyn_00', 'gsyn_01', 'gsyn_02', 'gsyn_03',
               'gsyn_10', 'gsyn_11', 'gsyn_12', 'gsyn_13',
               'gsyn_20', 'gsyn_21', 'gsyn_22', 'gsyn_23',
               'gsyn_30', 'gsyn_31', 'gsyn_32', 'gsyn_33',
               'stim0', 'stim1', 'stim2', 'stim3', 'stim4', 'stim5',
               'baseline', 'peak_thres', 'spike_0',
               'amp_0', 'slope_0', 'spike_1', 'amp_1', 'slope_1']

    # Dictionary of synapses categorized by population
    syn_pops = lib.get_synapses_gsyns_12()[1]

    # Get baseline of the inhibitory and excitatory neurons
    stims = [0*pamp]
    statemon = lib.run_model(eqs, reset_eqs, thres_eqs, syn_model, syn_eqs, syn_onpre,
                             params_e, params_i, variables, Vinit, stims, dt,
                             plot=False, deprived=deprived)[0]
    baseline = statemon[6][-1]
    baseline_i = statemon[9][-1]

    # Set Vinits to the baselines
    type_ind = {'exc': [0,1,2,6,7,8],
                'inh': [3,4,5,9,10,11]
    }
    Vinit = np.ones(N)
    for neuron in type_ind['exc']:
        Vinit[neuron] = baseline
    for neuron in type_ind['inh']:
        Vinit[neuron] = baseline_i

    # Set mean and sd for stimulation intensities
    stim_mu_e, stim_sd_e, stim_mu_i, stim_sd_i = 200*pamp, 15*pamp, 118*pamp, 15*pamp

    # Get N_iter results
    for times in range(N_iter):
        # Report progress
        if times%50 == 0:
            print 'Iteration', str(times+1) + '/' + str(N_iter)
        
        # Change gsyn in variables by random values per
        # pop->pop synapse, and get these values            
        variables['gsyn'], rand_gsyns = lib.get_random_gsyns(syn_pops, 0, max_gsyn, N)
        # Add gsyns to values_row
        values_row = []
        for i in range(4):
            for j in range(4):
                # Put value in row to write later
                values_row.append(rand_gsyns[(i,j)])         

        try:
            statemons = []
            # Simulate with different input currents avg_N times
            for avg_i in range(avg_N):
                # Change input currents
                stims = [np.random.normal(stim_mu_e, stim_sd_e) for i in range(N/4)] +\
                        [np.random.normal(stim_mu_i, stim_sd_i) for i in range(N/4)]

                # Run the network
                statemon, spikemon = lib.run_model(eqs, reset_eqs, thres_eqs, syn_model, syn_eqs,
                                        syn_onpre, params_e, params_i, variables, Vinit,
                                        stims, dt, N, runtime=100*ms, plot=plot, deprived=deprived)[0:2]
                statemons.append(statemon[6])

            # Get means of statemons
            statemons_zip = zip(*statemons)
            statemons_mean = np.asarray([mean(statemons_zip[i]) for i in range(len(statemons_zip))])

            # Get indices of the peaks
            peaks = lib.get_peaks(statemons_mean, baseline, peak_thres)
            # Get first amplitude
            amp = statemons_mean[peaks[0]] - baseline
            # Get first slope
            slope = lib.get_slope(statemon[6], peaks[0], baseline, .2)[0]

            # Save the second amp and slope too, because the first could just be a 'bump'
            if len(peaks) > 1:
                # Get amplitudes
                amp_1 = statemons_mean[peaks[1]] - baseline
                # Get slope
                slope_1 = lib.get_slope(statemon[6], peaks[1], baseline, .2)[0] 
            else:
                amp_1, slope_1 = None, None

            # Record whether the two peaks are spikes
            spike_0, spike_1 = False, False
            if len(spikemon[6]) > 0:
                if (len(peaks) > 0 and peaks[0] + spike_thres >= spikemon[6][0]/sampling_rate):
                    spike_0, spike_1 = True, True
                elif (len(peaks) > 1 and peaks[1] + spike_thres >= spikemon[6][0]/sampling_rate):
                    spike_0, spike_1 = False, True

            # Add results to values_row
            values_row += stims + [baseline, peak_thres, spike_0, amp, slope, spike_1, amp_1, slope_1]


            # If file does not yet exist, write columnheaders first 
            if os.path.isfile(output_file_name) is False:
                with open(output_file_name, "w") as f:
                    writer = csv.writer(f)
                    writer.writerow(headers)

            # Write values_row
            with open(output_file_name, "a") as f:
                writer = csv.writer(f)
                writer.writerow(values_row)


        # Exit program if a memoryerror occurs
        except MemoryError:
            print 'memory error'
            sys.exit()

        # If there is another error, save None as the amp/slope values
        except Exception as e:
            print e
            # Write the gsyns and None's
            values_row = values_row[:22] + [None, None, None, None]
            # If file does not yet exist, write columnheaders first
            if os.path.isfile(output_file_name) is False:
                with open(output_file_name, "w") as f:
                    writer = csv.writer(f)
                    writer.writerow(headers)

            # Write values_row
            with open(output_file_name, "a") as f:
                writer = csv.writer(f)
                writer.writerow(values_row)
