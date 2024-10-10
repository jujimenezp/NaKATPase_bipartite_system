#!/usr/bin/python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data = pd.read_csv("results/14_states_bw_reactions.dat", sep=";")
prop = data.iloc[:,0]
Edot_x = data.iloc[:,1]
Idot_x = data.iloc[:,2]
J = data.iloc[:,3]

fig, (ax1) = plt.subplots(1,1, layout='tight')

fig.suptitle("Energy and Information rates for subsystem X")
ax1.grid()
#ax2.grid()
ax1.set_xlabel(r"$\frac{W_{0,1}}{W_{7,8}}$")
#ax2.set_xlabel(r"$\frac{W_{0,1}}{W_{7,8}}$")
#ax1.set_ylabel(r"Energy rate ($k_\text{B}T/s$)")
#ax1.set_ylabel(r"Energy rate($k_\text{B}T/s$)")
#ax2.set_ylabel("Information rate (bits$/s$)")
ax1.set_ylim(None,0.017)
#ax1.set_xscale('log')
#ax2.set_xscale('log')
#ax2.set_xlim(0,500)
ax1.plot(prop,Edot_x, color='green', marker='.', linewidth=1, markersize=3, label=r'Energy rate ($k_\text{B}T/s$)')
ax1.plot(prop,Idot_x, color='red', marker='.', linewidth=1, markersize=3, label=r'Info. rate (bits$/s$)')
ax1.legend(loc='upper left')

plt.show()
#plt.savefig('results/bw_reactions_graph.png', format='png', bbox_inches='tight',dpi=300)
