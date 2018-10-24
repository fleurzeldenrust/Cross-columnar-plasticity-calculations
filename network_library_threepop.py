"""Library for network_threepop.py
with three neuron per population (12 total)

06-2018, Linda Wouters
"""

from brian import *
import pandas as pd
from scipy import stats


def get_synapses_gsyns_12():
    """Create a dictionary with synapse numbers
    categorized by connection type (i.e. excitatory
    to excitatory, exc. to inh., etc.)
    And a dictionary of synapses categorized by population"""
    N = 12
    # Neuron types
    ind_i_e = {'E': [0, 1, 2, 6, 7, 8],
               'I': [3, 4, 5, 9, 10, 11]}

    # Create dictionary of lists of synapses
    # sorted by type->type connection
    synapses = {"E_E": [],
                "E_I": [],
                "I_E": [],
                "I_I": []}
    # Iterate over each neuron
    for i in range(N):
        # If the neuron is excitatory, add tuples
        # to the E_x keys in synapses
        if i in ind_i_e['E']:
            synapses['E_E'] += [[i,x] for x in ind_i_e['E']]
            synapses['E_I'] += [[i,x] for x in ind_i_e['I']]
        # If the neuron is inhibitory, add tuples
        # to the I_x keys in synapses
        elif i in ind_i_e['I']:
            synapses['I_E'] += [[i,x] for x in ind_i_e['E']]
            synapses['I_I'] += [[i,x] for x in ind_i_e['I']]
        else:
            raise ValueError('Range neuronnumber outside\
                              of ind_i_e range') 


    # Create dictionary of lists of synapses
    # sorted by population->population connection:
    # syn_pops = {(0,0): [[0,0], [0,1], [0,2], ...],
    #              (0,1): [[0,3], [0,4], [0,5], ...], 
    #               ...}
    type_ind = {0: [0,1,2],
                1: [3,4,5],
                2: [6,7,8],
                3: [9,10,11]
    }
    syn_pops = {}
    # Iterate over all neuroncombinations and sort the synapses
    # according to population connection
    for key_pre in type_ind:
        for neur_pre in type_ind[key_pre]:

            for key_post in type_ind:
                for neur_post in type_ind[key_post]:
                    # Add neuron combi to correct connection type
                    if (key_pre, key_post) in syn_pops: 
                        syn_pops[(key_pre, key_post)].append([neur_pre, neur_post])
                    else:
                        syn_pops[(key_pre, key_post)] = [[neur_pre, neur_post]]

    return synapses, syn_pops


def get_random_gsyns(syn_pops, min, max, N):
    """Receives gsyns (dictionary of synapses categorized by populations)
    from get_synapses_gsyns_12() as input.
    Returns gsyns; a N*N array of gsyns, with each population->population
    synapse consisting of the same gsyn. 
    Returns gsyns, with gsyn for each synapse, 
    and rand_gsyns with the gsyn per connection type"""
    # get random gsyn values
    rand_gsyns = {key:np.random.uniform(min,max) for key in syn_pops}

    gsyns = np.ones((N,N))
    # Set gsyn values
    for key in syn_pops:
        for syn in syn_pops[key]:
            gsyns[syn[0]][syn[1]] = rand_gsyns[key]
    
    return gsyns, rand_gsyns


def get_random_Vinit(mu_e, sigma_e, mu_i, sigma_i, N):
    """Create list of initial membrane potentials,
    sorted by neuron index"""
    type_ind = {'exc': [0,1,2,6,7,8],
                'inh': [3,4,5,9,10,11],
    }
    Vinit = np.ones(N)

    # Set initial membrane potential for each neuron
    for neuron in type_ind['exc']:
        Vinit[neuron] = np.random.normal(mu_e, sigma_e)
    for neuron in type_ind['inh']:
        Vinit[neuron] = np.random.normal(mu_i, sigma_i)

    return Vinit


def get_initial_vals():
    """Returns the initial values for:
    - N (number of neurons)
    - eqs, reset_eqs, thres_eqs (neuron model equations)
    - syn_eqs, syn_onpre (synapse model equations)
    - syn_model (name of the synapse model)
    - params_e, params_i (paramter values for the exitatory
                            and inhibitory neuron models)
    - variables (dictionary of dictionaries, each describing
                the parameter values tau, Vs, and gsyn for
                all synapses
    - Vinit (list of initial membrane potential of the 12 neurons)
    - dt (1 / sampling rate)"""
    # Number of neurons
    N = 12
    # Get synapses categorized by type->type connections
    synapses = get_synapses_gsyns_12()[0]

    # Equations
    eqs = """
        dV/dt = ((-gL*(V-EL) + gL*dT*exp((V-Vt)/dT) - w + I - Isyn1) / C) : volt
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
        Isyn1 : amp
    """
    reset_eqs = """
        V = Vr
        w += b
    """
    thres_eqs = """
        V >= Vpeak
    """
    syn_eqs = """
    dx/dt = -x/tau_s : 1
    Isyn = gsyn * x * (V_post - Vs) : amp
    tau_s : second
    gsyn : siemens
    Vs : volt
    """
    syn_onpre = 'x+=1'

    # Parameters neuron model
    # Num 13, index 4
    params_e = {'a': 4.233569869733205e-10,
                'EL': -0.06602549741258842,
                'C': 9.392223810563386e-11,
                'b': 5.062387729563894e-12,
                'Vt': -0.04995878646849232,
                'tauw': 0.14283049763746267,
                'Vpeak': -0.00933816078033202,
                'dT': 0.0024175141183087976,
                'gL': 7.024557734819385e-09,
                'Vr': -0.062356325623533965
    } 
    # Num 1, index 12
    params_i = {'a': 3.541825376759559e-09,
                'EL': -0.06860671342824141,
                'C': 9.354222825336378e-11,
                'b': 3.76158704738044e-12,
                'Vt': -0.06167418618539455,
                'tauw': 0.1948899887860034,
                'Vpeak': 0.039695215712061103,
                'dT': 0.0010939061873941425,
                'gL': 1.4417940815623524e-08,
                'Vr': -0.08735205583835977
    }

    variables = {"tau_s": np.ones((N,N)),
                "Vs": np.ones((N,N)),
                "gsyn": np.ones((N,N))}
    # taus: http://homepages.inf.ed.ac.uk/mvanross/reprints/roth_mvr_chap.pdf
    # Vs: http://neuronaldynamics.epfl.ch/online/Ch3.S1.html
    type_vars = {
        "tau_s" : {"E_E": 1.7*ms,
                "E_I": .7*ms,
                "I_E": 6.5*ms,
                "I_I": 2.5*ms},
        "gsyn" : {"E_E": 10*nsiemens,
                "E_I": 10*nsiemens,
                "I_E": 10*nsiemens,
                "I_I": 10*nsiemens},
        "Vs" : {"E_E": 0*mV,
                "E_I": 0*mV,
                "I_E": -75*mV,
                "I_I": -75*mV}
    }
    # Create variables dict
    for var in type_vars:
        for key in type_vars[var]:
            for syn in synapses[key]:
                variables[var][syn[0]][syn[1]] = type_vars[var][key]

    # Initial membrane potentials
    Vinit = [-65*mV for __ in range(N)]

    syn_model = "exp"
    dt = .1*ms

    return eqs, reset_eqs, thres_eqs, syn_model, syn_eqs, syn_onpre, params_e, params_i, variables, Vinit, dt, N


def run_model(eqs, reset_eqs, thres_eqs, syn_model, syn_eqs,
              syn_onpre, params_e, params_i, variables, Vinit,
              stims, dt=.1*ms, N=12, runtime=100*ms, plot=False, deprived=False):
    """Uses values from get_initial_vals() as input:
    - eqs, reset_eqs, thres_eqs (neuron model equations)
    - syn_eqs, syn_onpre (synapse model equations)
    - syn_model (name of the synapse model)
    - params_e, params_i (paramter values for the exitatory
                            and inhibitory neuron models)
    - variables (dictionary of dictionaries, each describing
                the parameter values for E->E, E->I, I->E
                and I->I synapses (for tau, Vs, and gsyn))
    - Vinit (list of initial membrane potential of the four neurons)
    - stims (list of stimulus intensities for the neurons in the
             deprived column) 
    - dt (1 / sampling rate)
    - N (number of neurons)
    - runtime
    - plot (if True, model variables are plotted)
    - deprived (if True, Vt in excitatory neurons is changed)
    runs the fully connected model with these values and returns the spikemonitor and
    statemonitors (i.e. list of values per timestep, per neuron/synapse)
    for Vm, w, I, Isyn(neuron), Isyn(synapse) and xS"""

    defaultclock.reinit(t=0*second)
    defaultclock.dt = dt

    # Dictionary of neuron populations with lists
    # of the neuronindices corresponding to these
    # populations. With:
    # 0: excitatory deprived
    # 1: inhibitory deprived
    # 2: excitatory spared
    # 3: inhibitory spared
    type_ind = {0: [0,1,2],
                1: [3,4,5],
                2: [6,7,8],
                3: [9,10,11]
    }

    # Create network with synapses
    G = NeuronGroup(N, eqs, threshold=thres_eqs, reset=reset_eqs, clock=defaultclock)
    S = Synapses(G, G, model=syn_eqs,
                    pre=syn_onpre)
    # Variable name in the equations of G should have a different
    # name than the variable in S, or else it doesn't work
    G.Isyn1=S.Isyn
    S[:,:]=True

    # Set synaptic variables
    for i in range(N):
        for j in range(N):
            if syn_model == "exp":
                ind = S.synapse_index((i,j))
                S[ind[0]].tau_s = variables["tau_s"][i,j]
                S[ind[0]].Vs = variables["Vs"][i,j]
                S[ind[0]].gsyn = variables["gsyn"][i,j]

    # Initialize neuron variables and parameters
    for cat in [0, 2]:
        for i in type_ind[cat]:
            G[i].V = Vinit[i]
            G[i].C = params_e['C']
            G[i].gL = params_e['gL']
            G[i].EL = params_e['EL']
            G[i].Vt = params_e['Vt']
            G[i].dT = params_e['dT']
            G[i].a = params_e['a']
            G[i].b = params_e['b']
            G[i].tauw = params_e['tauw']
            G[i].Vr = params_e['Vr']
            G[i].Vpeak = params_e['Vpeak']
    for cat in [1, 3]:
        for i in type_ind[cat]:
            G[i].V = Vinit[i]
            G[i].C = params_i['C']
            G[i].gL = params_i['gL']
            G[i].EL = params_i['EL']
            G[i].Vt = params_i['Vt']
            G[i].dT = params_i['dT']
            G[i].a = params_i['a']
            G[i].b = params_i['b']
            G[i].tauw = params_i['tauw']
            G[i].Vr = params_i['Vr']
            G[i].Vpeak = params_i['Vpeak']

    # For deprived condition, change the threshold of spared exc 
    # (neuron population 2) with the same proportion as in the experimental
    # data (See excel sheet)
    if deprived is True:
        for i in type_ind[2]:
            G[i].Vt = params_e['Vt'] * (1-0.179598795)

    # Record spared pyramidal neuron (neuron_index=2)
    statemon = StateMonitor(G, 'V', record=True, clock=defaultclock)
    statemon_I = StateMonitor(G, 'I', record=True, clock=defaultclock)
    statemon_w = StateMonitor(G, 'w', record=True, clock=defaultclock)
    statemon_xS = StateMonitor(S, 'x', record=True, clock=defaultclock)
    statemon_Isyn = StateMonitor(S, 'Isyn', record=True, clock=defaultclock)
    statemon_Isyn1= StateMonitor(G, 'Isyn1', record=True, clock=defaultclock)
    spikemon = SpikeMonitor(G, record=True)

    net = Network(G, S, statemon, spikemon, statemon_w, statemon_I, statemon_Isyn, statemon_Isyn1, statemon_xS) 
 
    # Run model
    for i in range(len(stims)):
        G[i].I = stims[i]
    net.run(runtime)


    """Optional: plot the variables"""
    if plot == True:
        print "Spike times:"
        for i in range(N):
            print i, spikemon[i]

        plt.subplots_adjust(hspace=.8, wspace=.4)
        fig, ax = plt.subplots(3,2,figsize=(15,15))

        for i in range(N):
            ax[0,0].plot(statemon.times/ms, statemon[i]/mV, label=str(i))
        ax[0,0].legend()
        ax[0,0].set_title("V")

        for i in range(N):
            ax[0,1].plot(statemon_I.times/ms, statemon_I[i]/pamp, label=str(i))
        ax[0,1].legend()
        ax[0,1].set_title("I")

        for i in range(N):
            ax[1,0].plot(statemon_w.times/ms, statemon_w[i]/pamp, label=str(i))
        ax[1,0].legend()
        ax[1,0].set_title("w")

        for i in range(N*N):
            ax[1,1].plot(statemon_xS.times/ms, statemon_xS[i], label=str(i))
        ax[1,1].legend()
        ax[1,1].set_title("x in synapsemodel")

        for i in range(N*N):
            ax[2,0].plot(statemon_Isyn.times/ms, statemon_Isyn[i]/pamp, label=str(i))
        ax[2,0].legend()
        ax[2,0].set_title("Isyn in synapsemodel")

        for i in range(N):
            ax[2,1].plot(statemon_Isyn1.times/ms, statemon_Isyn1[i]/pamp, label=str(i))
        ax[2,1].legend()
        ax[2,1].set_title("Isyn in neuronmodel")

        plt.show()
    
    return statemon, spikemon, statemon_w, statemon_I, statemon_Isyn, statemon_Isyn1, statemon_xS


# def get_peaks_troughs(y):
#     try:
#         peaks = []
#         troughs = []
#         # Set the first difference
#         prev = y[1] - y[0]

#         # Iterate over each value and find the indices where the slope
#         # changes sign
#         for i in range(2, len(y)):
#             diff = y[i] - y[i-1]

#             # if previous difference was negative and current difference positive
#             # i-1 is a trough index
#             if prev < 0 and diff > 0:
#                 troughs.append(i-1)
#             # if previous difference was negative and current difference positive
#             # i-1 is a peak index
#             elif prev > 0 and diff < 0:
#                 peaks.append(i-1) 
#             # If prev is 0, no peak is stored, since we don't want the onset
#             # of the stimulation to count as a peak/trough

#             prev = diff

#         return peaks, troughs
#     except:
#         return None, None


def get_peaks(y, baseline, thres):
    """Get a list of peaks in y, with the peaks having
    an absolute amplitude (compared to baseline) bigger than thres
    If no peak is found, return None"""
    try:
        peaks = []
        # Set the first difference
        prev = y[1] - y[0]

        # Iterate over each value and find the indices where the slope
        # changes sign
        for i in range(2, len(y)):
            diff = y[i] - y[i-1]
            # if diff changes sign, return  i-1 as index of the peak
            if (((prev < 0 and diff > 0) or (prev > 0 and diff < 0))
                and abs(y[i-1]-baseline) > thres):
                peaks.append(i-1)
            # Set prev at current difference
            prev = diff

        # If prev is 0, no peak is stored, since we don't want the onset
        # of the stimulation to count as a peak/trough
        if len(peaks) > 0:
            return peaks
        else:
            return []
    except Exception as e:
        print e
        return None



def get_slope(y, peak_ind, baseline, prop, base_range=.00005):
    """ Takes as input the y value of the peak, the index of the peak,
    baseline, proportion (.2 in https://www.heka.com/downloads/software/manual/m_patchmaster.pdf p.178).
    Returns the slope, begin index and end index used for 
    calculating the slope"""
    try:
        # Determine the y values of the slope boundaries
        begin_y = baseline + (y[peak_ind] - baseline)*prop
        end_y = y[peak_ind] - (y[peak_ind] - baseline)*prop

        begin_ind = None
        end_ind = None
        # If peak is positive
        if y[peak_ind] > baseline:
            # Find the indices closest to begin and end y value
            for i in range(peak_ind, -1, -1):
                if end_ind is None and y[i] <= end_y:
                    end_ind = i
                    if begin_ind is not None:
                        break
                elif begin_ind is None and y[i] <= begin_y:
                    begin_ind = i
                    if end_ind is not None:
                        break
        # If peak is negative
        elif y[peak_ind] < baseline:
            # Find the indices closest to begin and end y value
            for i in range(peak_ind, -1, -1):
                # print i
                if end_ind is None and y[i] >= end_y:
                    end_ind = i
                    if begin_ind is not None:
                        break
                elif begin_ind is None and y[i] >= begin_y:
                    begin_ind = i
                    if end_ind is not None:
                        break
        # If the peak is equal to the baseline, return None
        else:
            return None, None, None
    
        # Return the slope of the linear regression between
        # begin_ind and end_ind
        return stats.linregress(range(begin_ind, end_ind+1), y[begin_ind:end_ind+1])[0], begin_ind, end_ind
    except Exception as e:
        print e
        return None, None, None
