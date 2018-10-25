"""Library for optimization with optim_multi.py
Optimized AdEx model to experimental data

Version including optimization of Vpeak variable
26-05-2018, Linda Wouters
"""

from brian import *
from brian.library.modelfitting import *
import pandas as pd
import time
import csv
import os
import neuronpy.util.spiketrain as spiketrain

def read_data(spike_file, current_file, sampling_rate, membr_file=None):
    """Takes as input a file of spikeindices and inputcurrent and
    the sampling rate.
    Returns numpy array of spike times and current input"""

    # Read spike indices data
    spikes = pd.read_csv(spike_file, header=None)
    current = pd.read_csv(current_file, header=None, names=["input"])
    spikes.transpose()

    # Divide spike indices by sampling rate (20) to get the spike times
    if membr_file is not None: 
        membr = pd.read_csv(membr_file, header=None, names=["Vm"])
        membr = membr["Vm"].values
    else:
        membr = None
    
    return spikes.iloc[0].values/sampling_rate, current["input"].values, membr



def write_results(model, results, output_file, input_file=None):
    """Write the dict results in a csv file"""
    
    # Set key_list and dimens for this model
    if model == "adex":
        # Sorted columnheaders
        key_list = ["best_fit", "C", "gL", "EL", "dT", "Vt", "a", "tauw", "b", "Vr", "Vpeak",
                    "spike_file", "current_length", "fit_time_ms", "data_dt",
                    "spikes_length", "popsize", "maxiter", "delta", "cpus",
                    "run_time", "date", "model", "neuron_type"]
        # Dimensions per parameter to remove units
        dimens = {"delta": 1*second, "C": 1*farad, "gL": 1*siemens, "EL": 1*volt, "dT": 1*volt, 
                  "Vt": 1*volt, "a": 1*siemens, "tauw": 1*second, "b": 1*amp, "Vr": 1*volt, "Vpeak": 1*volt}
    
    elif model == "izhik":
        # Sorted columnheaders        
        key_list = ["best_fit", "a", "b", "c", "d", "vt", "vr", "k", "C", "vpeak",
                    "spike_file", "current_length", "fit_time_ms", "data_dt",
                    "spikes_length", "popsize", "maxiter", "delta", "cpus",
                    "run_time", "date", "model", "neuron_type"]
        # Dimensions per parameter to remove units
        dimens = {"delta": 1*second, "a": 1/second, "b": 1*siemens, "c": 1*volt, "d": 1*amp, "vt": 1*volt, "vr": 1*volt, "k": 1*(siemens/volt), "C": 1*farad, "vpeak": 1*volt}
    else:
        raise ValueError("unknown model as input")

    # Create sorted row of the results
    values_row = []
    # Add None's for all values except spike_file
    if results is None:
        values_row = [input_file if key=="spike_file" else None for key in key_list]
    else:
        for key in key_list:
            # Remove dimension from values so it won't be written as a string
            if results[key] is not None and key in dimens:
                values_row.append(float(results[key]/dimens[key]))
            else:
                values_row.append(results[key])

    # If file does not yet exist, write columnheaders first
    if os.path.isfile(output_file) is False:
        with open(output_file, "w") as f:
            writer = csv.writer(f)
            writer.writerow(key_list)

    # Write values_row
    with open(output_file, "a") as f:
        writer = csv.writer(f)
        writer.writerow(values_row)


def modelfit(model, spikes, current, neuron_type, data_dt, popsize,
             maxiter, delta, cpus, spike_file, current_file, output_file):
    """Optimize data using Brians modelfitting function,
    returns a dict of results and writes the results using function
    'write_results()'"""
    assert type(spikes) is np.ndarray, "spikes is not a numpyarray: %r" % spikes
    assert type(current) is np.ndarray, "current is not a numpyarray: %r" % current
    assert (neuron_type == "exc" or neuron_type == "inh"), "neurontype is not exc or inh: %r" % neuron_type
    assert (data_dt > 0 and isinstance(data_dt, float))

    # Set dimensions
    spikes = spikes*ms
    current = current*pamp
    delta = delta*ms

    # Set dt
    defaultclock.dt = data_dt*ms

    # Get equations
    eqs, reset_eqs, thres_eqs = get_equations(model)
    # Get parameter ranges
    params = get_param_ranges(model, neuron_type)
    assert len(params) > 0

    # Optimize and time it
    set_time = time.time()
    results = modelfitting(model=eqs,
        reset=reset_eqs, threshold=thres_eqs,
        data=spikes, input=current, dt=defaultclock.dt, 
        popsize=popsize, maxiter=maxiter, delta=delta,
        cpu=cpus, **params)
    print results
    
    # Append best_fit and additional info to dict best_pos
    results.best_pos["best_fit"] = results.best_fit
    results.best_pos["popsize"] = popsize
    results.best_pos["maxiter"] = maxiter
    results.best_pos["delta"] = delta
    results.best_pos["cpus"] = cpus
    results.best_pos["spike_file"] = spike_file
    results.best_pos["current_file"] = current_file
    results.best_pos["spikes_length"] = len(spikes)
    results.best_pos["current_length"] = len(current)
    results.best_pos["fit_time_ms"] = len(current) * data_dt
    results.best_pos["data_dt"] = data_dt
    results.best_pos["model"] = model
    results.best_pos["neuron_type"] = neuron_type
    results.best_pos["run_time"] = time.time() - set_time
    results.best_pos["date"] = time.strftime("%x") + "_" + time.strftime("%X")

    print results.best_pos

    # Write results to output_file
    write_results(model, results, output_file)

    # Append results to output file
    with open("backup_results_dict.txt", "a+") as f:
        f.write(str(results.best_pos))
        # Seperate entries by newline
        f.write("\n")

    return results.best_pos


def get_equations(model):
    """Returns equations for 'adex' or 'izhik' model,
    takes model as input"""

    assert (model == "adex" or model == "izhik")

    if model == "adex":
        eqs = Equations("""
            dV/dt = ((-gL*(V-EL) + gL*dT*exp((V-Vt)/dT) - w + I) / C) : volt
            dw/dt = ((a*(V-EL) - w) / tauw) : amp
            I : amp
            C : farad
            gL : siemens
            EL : volt
            Vt : volt
            dT : volt
            a : siemens
            tauw : second
            b : amp
            Vr : volt
            Vpeak : volt
        """)
        reset_eqs = """
            V = Vr
            w += b
        """
        thres_eqs = """
            V >= Vpeak
        """
    
    elif model == "izhik":
        eqs = Equations("""
            dv/dt = (k*(v-vr)*(v-vt) - u + I) / C: volt
            du/dt = (a * (b*(v - vr) - u)): amp
            I : amp
            vt : volt
            a : 1 / second
            b : siemens
            c : volt
            d : amp
            vr : volt
            k : siemens/volt
            C : farad
            vpeak : volt
        """)
        reset_eqs = """
            v = c
            u = u + d
        """
        thres_eqs = """
            v >= vpeak
        """
    
    return eqs, reset_eqs, thres_eqs


def get_param_ranges(model, neuron_type):
    """Returns dictionary of parameter ranges
    for specific 'adex' or 'izhik' model
    and for 'exc' or 'inh' neuron type"""

    assert (model == "adex" or model == "izhik")
    assert (neuron_type == "exc" or neuron_type == "inh")

    if model == "adex" and neuron_type == "exc":
        # Broader parameter ranges based on min/max of
        # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2798047/#Sec4title
        # and the min/max values of previously optimized values
        # Same ranges as inhibitory adex
        params = {
            'C': [40*pF, 50*pF, 210*pF, 220*pF],
            'gL': [.05*nS, .1*nS, 28*nS, 30*nS],
            'EL' : [-90*mV, -85*mV, -35*mV, -20*mV],
            'dT' : [.1*mV, .5*mV, 10*mV, 20*mV],
            'Vt': [-75*mV, -70*mV, -20*mV, -10*mV],
            'a': [-12.5*nS, -11.75*nS, 6*nS, 7*nS],
            'tauw': [0*ms, 1*ms, 370*ms, 400*ms],
            'b': [0*pA, 0*pA, 130*pA, 160*pA],
            'Vr': [-95*mV, -90*mV, -50*mV, -40*mV],
            'Vpeak': [-20*mV, -10*mV, 40*mV, 50*mV]
        }
        return params
    
    elif model == "adex" and neuron_type == "inh":
        # Broader parameter ranges based on min/max of
        # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2798047/#Sec4title
        # and the min/max values of previously optimized values
        # Same ranges as excitatory adex
        params = {
            'C': [40*pF, 50*pF, 210*pF, 220*pF],
            'gL': [.05*nS, .1*nS, 28*nS, 30*nS],
            'EL' : [-90*mV, -85*mV, -35*mV, -20*mV],
            'dT' : [.1*mV, .5*mV, 10*mV, 20*mV],
            'Vt': [-75*mV, -70*mV, -20*mV, -10*mV],
            'a': [-12.5*nS, -11.75*nS, 6*nS, 7*nS],
            'tauw': [0*ms, 1*ms, 370*ms, 400*ms],
            'b': [0*pA, 0*pA, 130*pA, 160*pA],
            'Vr': [-95*mV, -90*mV, -50*mV, -40*mV],
            'Vpeak': [-20*mV, -10*mV, 40*mV, 50*mV]
        }
        return params
    
    elif model == "izhik" and neuron_type == "exc":
        # Based on min/max of article and short optimizations
        params = {
            'a': [0/ms, 0/ms, 0.25/ms, 0.3/ms],
            'c': [-100*mV, -90*mV, -40*mV, -20*mV],
            'b': [-1*nS, -0.7*nS, 0.35*nS, 0.45*nS],
            "vt": [-75*mV, -70*mV, 30*mV, 60*mV],
            'd': [-2.8*pA, -1.9*pA, 180*pA, 220*pA],
            'vr': [-90*mV, -80*mV, -50*mV, -40*mV],
            'k': [0*nS/mV, 0*nS/mV, .95*nS/mV, 1*nS/mV],
            'C': [1*pF, 10*pF, 150*pF, 200*pF],
            'vpeak': [-20*mV, -10*mV, 40*mV, 50*mV]
        }
        return params

    
    elif model == "izhik" and neuron_type == "inh":
        # Based on min/max of article and short optimizations     
        params = {
            'a': [0/ms, 0/ms, 0.25/ms, 0.3/ms],
            'c': [-100*mV, -90*mV, -40*mV, -20*mV],
            'b': [-1*nS, -0.7*nS, 0.35*nS, 0.45*nS],
            "vt": [-75*mV, -70*mV, 30*mV, 60*mV],
            'd': [-2.8*pA, -1.9*pA, 180*pA, 220*pA],
            'vr': [-90*mV, -80*mV, -50*mV, -40*mV],
            'k': [0*nS/mV, 0*nS/mV, .95*nS/mV, 1*nS/mV],
            'C': [1*pF, 10*pF, 150*pF, 200*pF],
            'vpeak': [-20*mV, -10*mV, 40*mV, 50*mV]
        }
        return params


def run_model(model, result, current, SAMPLE_RATE, model_time, Vpeak=None):
    """Run the model with results as parameters.
    Returns the spiketimes, membrane potential and modeled times"""
    assert (model == "adex" or model == "izhik")
    
    # Set defaultclock
    defaultclock.reinit(t=0*second)
    print SAMPLE_RATE
    print float(SAMPLE_RATE)
    defaultclock.dt = (1/float(SAMPLE_RATE)) * ms

    # Set equations, network and parameters
    if model == "adex":
        eqs, reset_eqs, thres_eqs = get_equations("adex")
        # params = list(get_param_ranges(model, neuron_type))
        # Network
        G = NeuronGroup(1, eqs, threshold=thres_eqs, reset=reset_eqs, clock=defaultclock)

        # Give input
        G.I = TimedArray(current * pamp, dt=defaultclock.dt)
        # Set parameters
        G.C = result['C']
        G.gL = result['gL']
        G.EL = result['EL']
        G.Vt = result['Vt']
        G.dT = result['dT']
        G.a = result['a']
        G.b = result['b']
        G.tauw = result['tauw']
        G.Vr = result['Vr']
        G.Vpeak = result['Vpeak']

        # Monitors for amplitude and spikes
        statemon = StateMonitor(G, 'V', record=0, clock=defaultclock)
        spikemon = SpikeMonitor(G)


    elif model == "izhik":
        eqs, reset_eqs, thres_eqs = get_equations("izhik")

        # Network
        G = NeuronGroup(1, eqs, threshold=thres_eqs, reset=reset_eqs, clock=defaultclock)

        # Give input
        G.I = TimedArray(current * pamp, dt=defaultclock.dt)

        # Set parameters
        G.a = result['a']
        G.b = result['b']
        G.c = result['c']
        G.d = result['d']
        G.vt = result['vt']
        G.vr = result['vr']
        G.k = result['k']
        G.C = result['C']
        G.vpeak = result['vpeak']

        # Monitors for amplitude and spikes
        statemon = StateMonitor(G, 'v', record=0, clock=defaultclock)
        spikemon = SpikeMonitor(G)

    run(model_time*ms)

    # Return spiketimes and Vm
    return spikemon.spiketimes[0]/ms, statemon[0]/mV, statemon.times/ms 


def differences(modeled_vm, experimental_vm):
    """Returns the MSE and difference between the two means"""

    # Set length
    if len(modeled_vm) > len(experimental_vm):
        length = len(experimental_vm)
    else:
        length = len(modeled_vm)

    # Calculate mean square estimate
    MSE = sum([(modeled_vm[i] - experimental_vm[i])**2 for i in range(length)]) / length

    return MSE, mean(modeled_vm), mean(experimental_vm)


def write_gammas(modeled_spikes, spikes, result, model_time, spike_file, output_file, model, modeled_vm, experimental_vm):
    """Write results of the retest
    Takes as input:
    - modeled spikes
    - spikes: experimental spikes
    - result: the first-optimization result
    - model_time: time the model is run
    - spike_file: the name of the file containing the experimental files
    - output_file: name of the output file
    - model: modeltype ('izhik' or 'adex')
    - modeled_vm: modeled membrane potential
    - experimental_vm: expermental membrane potential.
    Returns the row of values that are written into the output file."""

    # Get the experimental and modeled spikes within the 'retest' timeframe
    retest_data_spikes = [x for x in spikes if (x >= result["fit_time_ms"] and x < model_time)]
    retest_modeled_spikes = [x for x in modeled_spikes if (x >= result["fit_time_ms"] and x < model_time)]

    # Set column headers and first values
    # Retest time: the time corresponding to the amount of data 
    column_names = ["retest_time", "retest_data_spike_length", "retest_modeled_spike_length",
                    "spike_file", "optim_fit", "model", "vm_MSE", "modeled_vm_mean",
                    "exp_vm_mean"]
    values_row = [model_time - result["fit_time_ms"], len(retest_data_spikes),
                  len(retest_modeled_spikes), spike_file, result["best_fit"],
                  model]

    # Calculate the MSE and mean modeled and experimental membrane potential
    diff = differences(modeled_vm, experimental_vm)
    values_row += [diff[0], diff[1], diff[2]]

    # Get coincidence factor (gamma) and number of coincidences
    for i in range(1, 11):
        column_names += ["gamma_%s" %str(i), "coincidences_%s" %str(i)]
        values_row.append(coincidence_factor(retest_data_spikes, retest_modeled_spikes, window=i))
        values_row.append(spiketrain.coincident_spikes(retest_data_spikes, retest_modeled_spikes, window=i, normalize=True)[1])


    # If file does not yet exist, write columnheaders first
    if os.path.isfile(output_file) is False:
        with open(output_file, "w") as f:
            writer = csv.writer(f)
            writer.writerow(column_names)

    # Write values_row
    with open(output_file, "a") as f:
        writer = csv.writer(f)
        writer.writerow(values_row)

    return values_row


def coincidence_factor(ref, comp, window=5, isi=None):
    """ 
    https://pythonhosted.org/neuronpy/_modules/neuronpy/util/spiketrain.html#coincidence_factor
    coincidence factor from neuronpy's util.spiketrain with added "*(1 / (1 - 2 * v * window))"
    to get the same calculation of gamma as in Jolivet et al (2008)
    (http://icwww.epfl.ch/~gerstner/PUBLICATIONS/Jolivet08.pdf)

    The coincidence factor :math:`\Gamma` between two spike trains is defined as

    .. math::
        
       \Gamma = \frac{N_\mathrm{coinc}- E \left( N_\mathrm{coinc} \right)}
       {\frac{1}{2}\left(N_\mathrm{ref}+N_\mathrm{comp}\right) - 
       E \left( N_\mathrm{coinc} \right)}

    where :math:`N_{\mathrm{ref}}` are the number of spikes in the reference train,
    :math:`N_{\mathrm{comp}}` is the number of spikes in the comparing train, 
    :math:`N_{\mathrm{coinc}}` is the number of coincident spikes within a time window 
    :math:`\Delta`, :math:`E \left( N_\mathrm{coinc} \right) = 2 v \Delta N_{\mathrm{ref}}` 
    is the expected number of coincident spikes that would be given by chance 
    if the spikes in the comparing train were generated by a homogeneous 
    Poisson process with its rate :math:`v`. This correlation measure has the range 
    [-1, 1] where 1 is perfectly correlated, 0 is not correlated, and -1 is 
    perfectly anti-correlated.
    
    :param ref: Spike times of the reference train. It is 
        assumed that the spike times are ordered sequentially.
    
    :param comp: Spike times of the train that shifts in time.
        It is assumed that the spike times are ordered sequentially.
    
    :param window: Time window to say that two spikes are synchronized.
        This has a default value of 5.
        
    :param isi: If supplied, this is the isi of the comparing train. Otherwise,
        the rate of the train is computed by taking the last spike minus the
        first spike and dividing by the number of spikes in between.
    
    :return: Coincidence factor
    """
    num_coincident, mask_a, mask_b, ratio = spiketrain.get_sync_traits( \
            ref, comp, window)
    len_ref = len(ref)
    len_comp = len(comp)
    total_spikes = len_ref + len_comp
    if isi is None:
        v = (len_comp - 1)/(comp[-1] - comp[0])
    else:
        v = 1./isi
    expected_coincidences = 2 * v * window * len_ref
    # Added "*(1 / (1 - 2 * v * window))"
    return (num_coincident - expected_coincidences)*2/ \
            (total_spikes - (2*expected_coincidences))* \
            (1 / (1 - 2 * v * window))


def initialize_data(spike_file, min_spikes, SAMPLE_RATE):
    """Initialize the variables"""
    # Get file names from spike_file
    current_file = spike_file.replace("spikeindices", "input")
    membr_file = spike_file.replace("spikeindices", "membrane")

    # Get spike times and current input
    spikes, current, membrane = read_data(spike_file,
                                    current_file,
                                    int(SAMPLE_RATE), membr_file)
    # Downsample to half the sample rate
    current = current[::2]
    # Set downsampled sample_rate
    data_dt = (2/float(SAMPLE_RATE))

    # Check if there are enough spikes to optimize
    if len(spikes) < min_spikes:
        raise ValueError("Not enough spikes in file.")
    # Set range of optimization
    elif len(spikes)*2/3 > 51:
        # Set fitrange to 50th spike
        fit_range = spikes[50] / data_dt
    else:
        # Set fitrange to the 2/3rd spikelength-spike
        fit_range = spikes[int(len(spikes)*2/3)] / data_dt
    
    # Shorten data to speed up optimization
    spikes = np.array([x for x in spikes if x < fit_range*data_dt])
    current = current[:int(fit_range)]

    return spikes, current, membrane, data_dt, fit_range
