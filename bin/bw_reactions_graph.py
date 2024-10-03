#!/usr/bin/python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data = pd.read_csv("results/14_states_bw_reactions.dat", sep=";")
prop = data.iloc[:,0]
Edot_x = data.iloc[:,1]
Idot_x = data.iloc[:,2]

fig, (ax1, ax2) = plt.subplots(1,2, layout='tight')

ax1.grid()
ax2.grid()
fig.suptitle("Energy and Information rates for subsystem X")
ax1.set_xlabel(r"$\frac{W_{0,1}}{W_{7,8}}$")
ax2.set_xlabel(r"$\frac{W_{0,1}}{W_{7,8}}$")
ax1.plot(prop,Edot_x, color="green", label="Energy rate (X)")
ax2.plot(prop,Idot_x, color="red", label="Information rate (X)")
ax1.legend()
ax2.legend()

#plt.show()
plt.savefig('results/bw_reactions_graph.png', format='png', bbox_inches='tight', dpi=300)
