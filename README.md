# Cross-columnar-plasticity-calculations
Code for the fitting of neural data, network simuilations and analysis for the cross-columnar plasticity upon whisker deprivation described in (TOADD)

## Anaconda virtual environment
We provide an Anaconda virtual environment with all the packages that should make the installation of the right (older) versions of python (2.7.15) and numpy (1.10.4). For the modelfitting, playdoh (https://github.com/rossant/playdoh) is needed. For data storage and analysis, we used pandas (https://pandas.pydata.org/).

## Modelfitting
An Adaptive Exponential Integrate-and-fire neuron (see http://www.scholarpedia.org/article/Adaptive_exponential_integrate-and-fire_model) was fitted to in-vitro data (see TOADD), using the Brian 1 modelfitting toolbox (see https://brian.readthedocs.io/en/1.4.3/modelfitting.html, it has not been exported to Brian 2 yet). We provide an anaconda virtual environment to get the right installation. The python file optim_vpeak.py fits a AdEx model to all data files in the current directory. It needs library optimize_vpeak.py. To export matlab data to .csv files, the MATLAB file 'export_data.m' is provided.

## Network simulations
The network simulations can be done by entering either 'python network_threepop_Fleur.py control' or 'python network_threepop_Fleur.py deprived' into the terminal of the virtual environment.  The library 'network_library_threepop.py' is needed for these simulations. 

## Network selection
The selection of the networks was done according to the steps described in the Jupyter notebook TOADD
