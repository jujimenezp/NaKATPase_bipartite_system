# Analisys of the Na+,K+-ATPase

I performed an analysis on the 19x19 matrix constructed using the equations from Clarke et al. (2013) for the Albers-Post model that describes the behavior of the Na+,K+-ATPase, the protein responsible for the Sodium Potassium active transport accross the cellular membrane.
The explicit W matrix can be seen [here](data/Sodium_Potassium_pump_W_matrix.pdf).

A bipartite assumption was made for a simplified model considering only the main path of the protein. Rate of change of heat, work, and energy were calculated for the subsystems and the whole system as well.

This program can be run with the command
`python bin/driver.py`
and it requires the Eigen library for C++.

The core of the program is in the `src/` folder, namely, in the files W_matrix.hpp and 14_states.cpp. The file 19_states.cpp is outdated and is only kept for reference purposes.

The `data/` folder contains initial data such as transition rates and parameters for the system. Results are kept in the `results/` folder.

Additionally, in `bin/` there exist the programs sym_19_states.jl and sym_19_states.py that were an attempt to use symbolic algebra to find an analytical expression for the nonequilibrium steady state distribution. However, this is not operational at the moment.

R. J. Clarke et al., Biochimica et Biophysica Acta (BBA) - Bioenergetics 1827, 1205 (2013). 
