#!/usr/bin/python3

from pandas import read_csv
from matplotlib import pyplot as plt

data = read_csv("results/proportion_range.dat", sep="\t")

prop = data.iloc[:,0]
turnover = data.iloc[:,1]
Idot_x = data.iloc[:,8]
Efficiency = data.iloc[:,10]
Qdot_x = data.iloc[:,4]
Qdot_y = data.iloc[:,5]
Wdot_x = data.iloc[:,6]
Wdot_y = data.iloc[:,7]
Work = data.iloc[:,2]
Heat = data.iloc[:,3]
Edot_x = Wdot_x + Qdot_x

fig, ax1 = plt.subplots(1,1, figsize=(10,5), layout='tight')
#ax.set_title(r"Proportion v. Information Flow")
ax1.set_ylim(None,17)
ax1.set_xscale('log')
ax1.grid()
ax1.set_xlabel(r"$\frac{W_{0,1}}{W_{7,8}}$")
ax1.set_ylabel(r'$\Delta E_X$ ($k_\text{B}T/$cycle)', color='green')
ax1.plot(prop, Edot_x, color='green', linewidth=1)
ax1.tick_params(axis='y', labelcolor='green')
ax2 = ax1.twinx()
ax2.set_ylabel(r'$\Delta I_X$ (nats/cycle)', color='red')
ax2.plot(prop,Idot_x, color='red', linewidth=1)
ax2.tick_params(axis='y', labelcolor='red')
#plt.show()
#plt.savefig('results/prop_energy_info.png', format='png', bbox_inches='tight', dpi=300)

fig, ax = plt.subplots(1,1, layout='tight')
ax.set_xscale('log')
ax.grid()
ax.set_xlabel(r"$\frac{W_{0,1}}{W_{7,8}}$")
ax.set_ylabel('Probability current (1/s)')
ax.plot(prop, turnover, color='black', linewidth=1)
plt.show()
#plt.savefig('results/prop_turnover.png', format='png', bbox_inches='tight', dpi=300)
