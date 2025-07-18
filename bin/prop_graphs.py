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
eff_x = (Edot_x+Idot_x)/Wdot_x
eff_y = -Wdot_y/(Edot_x+Idot_x)

fig, (ax1, ax3) = plt.subplots(2,1,  sharex=True)
#ax.set_title(r"Proportion v. Information Flow")

#ax1.ticklabel_format(axis='y', style='scientific', scilimits=(-1,-1), useMathText=True)
#ax1.yaxis.set_major_locator(plt.MultipleLocator(2e-1))
ax1.set_yticks([0, 0.2, 0.4, 0.6, 0.8], ['0', '0.2', '0.4', '0.6', '0.8'])
ax1.set_xlim(1e-9,1e+9)
ax1.set_ylim(0,0.8)
ax1.set_box_aspect(0.6)
ax1.set_xscale('log')
#ax1.set_xlabel(r"$\frac{W_{0,1}}{W_{7,8}}$")
ax1.set_ylabel('Probability current (1/s)', fontsize=11)
ax1.plot(prop, turnover, color='black', linewidth=2)
ax1.text(0.05,0.95, 'a.', transform=ax1.transAxes, fontsize=13, verticalalignment='top', backgroundcolor='None')

#ax2.set_box_aspect(0.6)
#ax2.set_xlim(1e-9,1e+9)
#ax2.set_ylim(0,1)
#ax2.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1], ['0', '0.2', '0.4', '0.6', '0.8', '1'])
#ax2.set_xscale('log')
#ax2.set_ylabel(r'Efficiency', fontsize=11)
#ax2.plot(prop, Efficiency, color='purple', linewidth=2)
#ax2.plot(prop, eff_x, color='red', linewidth=2)
#ax2.plot(prop, eff_y, color='blue', linewidth=2)
#ax2.text(0.05,0.15, 'b.', transform=ax2.transAxes, fontsize=13, verticalalignment='top', backgroundcolor='None')
#ax2.text(0.8,0.85, r'$\eta_Y$', transform=ax2.transAxes, fontsize=11, verticalalignment='top', backgroundcolor='None', color='blue')
#ax2.text(0.8,0.99, r'$\eta_X$', transform=ax2.transAxes, fontsize=11, verticalalignment='top', backgroundcolor='None', color='red')
#ax2.text(0.45,0.75, r'$\eta$', transform=ax2.transAxes, fontsize=11, verticalalignment='top', backgroundcolor='None', color='purple')

ax3.set_box_aspect(0.6)
ax3.set_xlim(1e-9,1e+9)
ax3.set_ylim(None,20)
ax3.set_xscale('log')
ax3.set_xlabel(r"$\gamma$", fontsize=11)
ax3.set_ylabel('Thermodynamic flows (/cycle)', fontsize=10)
ax3.plot(prop, Edot_x, color='green', linewidth=2)
ax3.plot(prop,Idot_x, color='red', linewidth=2)
ax3.text(0.05,0.95, 'c.', transform=ax3.transAxes, fontsize=14, verticalalignment='top', backgroundcolor='None')
E_label = r'Energy $\beta \Delta E_X$'
I_label = r'Information $\Delta I_X$'
ax3.text(0.55,0.76, E_label, transform=ax3.transAxes, fontsize=8, verticalalignment='top', backgroundcolor='None', color='green')
ax3.text(0.55,0.38, I_label, transform=ax3.transAxes, fontsize=8, verticalalignment='top', backgroundcolor='None', color='red')


plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
extent = mpl.transforms.Bbox([[0.2,0.02],[0.8,0.472]]).transformed(fig.transFigure - fig.dpi_scale_trans)
plt.savefig('results/prop_energy_info.png', format='png', bbox_inches=extent, dpi=300)