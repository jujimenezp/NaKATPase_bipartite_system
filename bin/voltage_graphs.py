#!/usr/bin/python3

from pandas import read_csv
from matplotlib import pyplot as plt
import matplotlib as mpl
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'

data = read_csv("results/voltage_range.dat", sep="\t")

voltage = data.iloc[:,0]*1000
turnover = data.iloc[:,1]
Idot_x = data.iloc[:,8]
Efficiency = data.iloc[:,10]
Qdot_x = data.iloc[:,4]
Qdot_y = data.iloc[:,5]
Wdot_x = data.iloc[:,6]
Wdot_y = data.iloc[:,7]
Work = data.iloc[:,2]
Heat = data.iloc[:,3]
Energy = Work + Heat
Edot_x = Wdot_x + Qdot_x
eff_x = (Edot_x + Idot_x)/Wdot_x
eff_y = -Wdot_y/(Edot_x + Idot_x)

fig, ax = plt.subplots(1,1, layout='tight')

#ax.set_title("Voltage v. Turnover rate")
ax.grid()
ax.set_yscale('log')
ax.set_xlabel('Voltage (V)')
ax.set_ylabel('Turnover rate (1/s)')
ax.plot(voltage,turnover, color='black', marker='.', linewidth=1, markersize=3)
plt.savefig('results/voltage_turnover.png', format='png', bbox_inches='tight', dpi=300)

fig, ax = plt.subplots(1,1, layout='tight')
#ax.set_title(r"Voltage v. Information Flow")
ax.grid()
ax.set_xlabel('Voltage (V)')
ax.set_ylabel(r'$\Delta I_X$ (nats/cycle)')
ax.axhline(0, color='black', linestyle='dashed')
ax.plot(voltage,Idot_x, color='red', linewidth=1)
plt.savefig('results/voltage_info_flow.png', format='png', bbox_inches='tight', dpi=300)

fig, ax = plt.subplots(1,1, layout='tight')
ax.set_title("Voltage v. Efficiency")
ax.grid()
ax.set_xlabel('Voltage (V)')
ax.set_ylabel('Efficiency')
ax.plot(voltage,Efficiency, color='black', marker='.', linewidth=1, markersize=3)
plt.savefig('results/voltage_efficiency.png', format='png', bbox_inches='tight', dpi=300)

fig, ax1 = plt.subplots(1,1, layout='tight')
#ax1.set_title("Voltage v. Subsystems Heat")
ax1.grid()
ax1.set_xlabel('Voltage (V)')
ax1.set_ylabel(r'Heat ($\frac{k_B T}{cycle}$)')
ax1.axhline(0, color='black', linestyle='dashed')
ax1.plot(voltage,Qdot_x, color='red', linewidth=1, label=r'$\Delta Q_X$')
ax1.plot(voltage,Qdot_y, color='blue', linewidth=1, label=r'$\Delta Q_Y$')
ax1.legend(loc="upper right")
plt.savefig('results/voltage_sub_heat.png', format='png', bbox_inches='tight', dpi=300)

fig, ax2 = plt.subplots(1,1, layout='tight')
#ax2.set_title("Voltage v. Subsystems Work")
ax2.grid()
ax2.set_xlabel('Voltage (V)')
ax2.set_ylabel(r'Work ($\frac{k_B T}{cycle}$)')
ax2.plot(voltage,Wdot_x, color='red', linewidth=1, label=r'$\dot{W}_X$')
ax2.plot(voltage,Wdot_y, color='blue', linewidth=1, label=r'$\dot{W}_Y$')
ax2.legend(loc="upper right")
plt.savefig('results/voltage_sub_work.png', format='png', bbox_inches='tight', dpi=300)


fig, ax = plt.subplots(1,1, layout='tight')
ax.set_title("Voltage v. Heat and Work")
ax.grid()
ax.set_xlabel('Voltage (V)')
ax.set_ylabel(r'($\frac{k_B T}{cycle}$)')
ax.plot(voltage,Heat, color='black', marker='.', linewidth=1, markersize=3, label='Heat')
ax.plot(voltage,Work, color='red', marker='.', linewidth=1, markersize=3, label='Work')
ax.legend(loc='upper left')
plt.savefig('results/voltage_global.png', format='png', bbox_inches='tight', dpi=300)
#plt.show()


fig, (ax1, ax2) = plt.subplots(2,1, sharex=True)
ax1.text(0.05, 0.95, 'a.', transform=ax1.transAxes, fontsize=14, verticalalignment='top', backgroundcolor='None')
ax1.set_xlim(-0.087,0.06)
ax1.set_ylim(0,5.6e-3)
ax1.ticklabel_format(axis='y', style='scientific', scilimits=(-3,-3), useMathText=True)
ax1.set_box_aspect(0.6)
ax1.set_ylabel('Probability current (1/s)')
ax1.plot(voltage,turnover, color='black', linewidth=2)

# #ax1.set_xlim(-0.087,0.06)
# ax1.text(0.05, 0.95, 'b.', transform=ax1.transAxes, fontsize=14, verticalalignment='top', backgroundcolor='None')
# ax1.set_ylim(0,1)
# ax1.set_box_aspect(0.6)
# ax1.set_ylabel(r'Efficiency', fontsize=11)
# ax1.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1], ['0', '0.2', '0.4', '0.6', '0.8', '1'])
# ax1.plot(voltage, Efficiency, color='purple', linewidth=2, label=r'$\eta$')
# ax1.plot(voltage, eff_x, color='red', linestyle='dashed', linewidth=2, label=r'$\eta_X$')
# ax1.plot(voltage, eff_y, color='blue', linewidth=2, label=r'$\eta_Y$')
# #ax1.legend(loc="lower left")
# ax1.text(0.4, 0.98, r'$\eta_Y$', fontsize=14, transform_rotates_text=True, transform=ax1.transAxes, verticalalignment='top', backgroundcolor='None', color='blue')
# ax1.text(0.4, 0.82, r'$\eta_X$', fontsize=14, transform_rotates_text=True, transform=ax1.transAxes, verticalalignment='top', backgroundcolor='None', color='red')
# ax1.text(0.7, 0.6, r'$\eta$', fontsize=14, transform_rotates_text=True, transform=ax1.transAxes, verticalalignment='top', backgroundcolor='None', color='purple')

ax2.set_xlim(-87,60)
ax2.set_box_aspect(0.6)
ax2.set_xlabel('Voltage (mV)')
ax2.set_ylabel('Heat and information flow')
ax2.axhline(0, color='black', linestyle='dashed', linewidth=0.7)
ax2.plot(voltage,Qdot_x, color='red', linewidth=2, label=r'$\Delta Q_X$')
ax2.plot(voltage,Idot_x, color='orange', linewidth=2, label=r'$\Delta I_X$')
ax2.plot(voltage,Qdot_y, color='blue', linestyle='dashed', linewidth=1.5, label=r'$\Delta Q_Y$')
ax2.plot(voltage, Edot_x, color='green', linewidth=2, label=r'$\Delta E_X$')
#ax2.legend(loc='best', fontsize=8, framealpha=1)

ax2.text(0.05, 0.95, 'c.', transform=ax2.transAxes, fontsize=14, verticalalignment='top', backgroundcolor='None')
ax2.text(0.1, 0.16, r'$\Delta Q_X$', rotation=13, rotation_mode='anchor', transform_rotates_text=True, transform=ax2.transAxes, verticalalignment='top', backgroundcolor='None', color='red')
ax2.text(0.1, 0.56, r'$\Delta I_X$', rotation=-16, rotation_mode='anchor', transform_rotates_text=True, transform=ax2.transAxes, verticalalignment='top', backgroundcolor='None', color='orange')
ax2.text(0.42, 0.32, r'$\Delta Q_Y$', rotation=-18, rotation_mode='anchor', transform_rotates_text=True, transform=ax2.transAxes, verticalalignment='top', backgroundcolor='None', color='blue')
ax2.text(0.35, 0.94, r'$\Delta E_X$', rotation=13, rotation_mode='anchor', transform_rotates_text=True, transform=ax2.transAxes, verticalalignment='top', backgroundcolor='None', color='green')


#fig.subplots_adjust(hspace=.0)
plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)

extent = mpl.transforms.Bbox([[0.2,0.515],[0.8,0.92]]).transformed(fig.transFigure - fig.dpi_scale_trans)
fig.savefig('results/voltage_turnover.png', format='png', bbox_inches=extent, dpi=300)
