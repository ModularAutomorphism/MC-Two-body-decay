# Simulating Two Body Decay process with two decay channels 
## Decay considered: Charged muon decay to muon + muon-antineutrino and electron (or positron) + electron-antineutrino

### C++ file for decay simulation
The function charged_pion_decay_simulator(int N) takes input the number of pions "N" to decay. By default the code file has N = 1000000. \
**Note:** Please simulate at least 10000 pions to prevent segmentation fault. This is because probability of decay to electron channel is very small and hence simulating less number of pions can result in empty arrays for electron channel data which will result in the faults.

### Jupyter Notebook for Visualization
The simulated data can be visualized using the Python code in the Jupyter notebook.
