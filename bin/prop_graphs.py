#!/usr/bin/python3

from pandas import read_csv
from matplotlib import pyplot as plt
import matplotlib as mpl
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'

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

fig, (ax1, ax2) = plt.subplots(2,1,  sharex=True)
#ax.set_title(r"Proportion v. Information Flow")

#fig, ax = plt.subplots(1,1, layout='tight')
#ax1.ticklabel_format(axis='y', style='scientific', scilimits=(0,0))
ax1.set_xlim(1e-9,1e+9)
ax1.set_box_aspect(0.6)
ax1.set_xscale('log')
#ax1.grid()
#ax1.set_xlabel(r"$\frac{W_{0,1}}{W_{7,8}}$")
ax1.set_ylabel('Probability current (1/s)', fontsize=11)
ax1.plot(prop, turnover, color='black', linewidth=2)
#plt.savefig('results/prop_turnover.png', format='png', bbox_inches='tight', dpi=300)

ax2.set_box_aspect(0.6)
ax2.set_xlim(1e-9,1e+9)
ax2.set_ylim(None,20)
ax2.set_xscale('log')
#ax2.grid()
ax2.set_xlabel(r"$\frac{R_{0,1}}{R_{7,8}}$", fontsize=13)
ax2.set_ylabel('Energy and Info. rates', fontsize=11)
ax2.plot(prop, Edot_x, color='green', linewidth=2, label=r'Energy rate ($k_\text{B}T/s$)')
ax2.plot(prop,Idot_x, color='red', linewidth=2,label=r'Info. rate (nats$/s$)')
ax2.legend(loc="upper right")
ax2.yaxis.get_major_ticks()[-1].label1.set_visible(False)

fig.subplots_adjust(hspace=.0)
plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
plt.savefig('results/prop_graph.png', format='png', bbox_inches='tight',dpi=300)
#plt.show()
#ax1.tick_params(axis='y', labelcolor='green')
#ax2 = ax1.twinx()
#ax2.set_ylabel(r'$\Delta I_X$ (nats/cycle)', color='red')
#ax2.plot(prop,Idot_x, color='red', linewidth=1)
#ax2.tick_params(axis='y', labelcolor='red')
#plt.savefig('results/prop_energy_info.png', format='png', bbox_inches='tight', dpi=300)
