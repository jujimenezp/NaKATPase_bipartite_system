#!/usr/bin/python3

from pandas import read_csv
from matplotlib import pyplot as plt

data = read_csv("results/voltage_range.dat", sep="\t")

voltage = data.iloc[:,0]
turnover = data.iloc[:,1]
Idot_x = data.iloc[:,8]
Efficiency = data.iloc[:,10]
Qdot_x = data.iloc[:,4]
Qdot_y = data.iloc[:,5]
Wdot_x = data.iloc[:,6]
Wdot_y = data.iloc[:,7]
Work = data.iloc[:,2]
Heat = data.iloc[:,3]

fig, ax = plt.subplots(1,1, layout='tight')

ax.set_title("Voltage v. Turnover rate")
ax.grid()
ax.set_yscale('log')
ax.set_xlabel('Voltage (V)')
ax.set_ylabel('Turnover rate (1/s)')
ax.plot(voltage,turnover, color='black', marker='.', linewidth=1, markersize=3)
plt.savefig('results/voltage_turnover.png', format='png', bbox_inches='tight', dpi=300)

fig, ax = plt.subplots(1,1, layout='tight')
ax.set_title(r"Voltage v. Information Flow")
ax.grid()
ax.set_xlabel('Voltage (V)')
ax.set_ylabel(r'$\dot{I}_X$ (nats/cycle)')
ax.plot(voltage,Idot_x, color='black', marker='.', linewidth=1, markersize=3)
plt.savefig('results/voltage_info_flow.png', format='png', bbox_inches='tight', dpi=300)

fig, ax = plt.subplots(1,1, layout='tight')
ax.set_title("Voltage v. Efficiency")
ax.grid()
ax.set_xlabel('Voltage (V)')
ax.set_ylabel('Efficiency')
ax.plot(voltage,Efficiency, color='black', marker='.', linewidth=1, markersize=3)
plt.savefig('results/voltage_efficiency.png', format='png', bbox_inches='tight', dpi=300)

fig, (ax1, ax2) = plt.subplots(1,2, layout='tight')
ax1.set_title("Voltage v. Subsystems Heat")
ax1.grid()
ax1.set_xlabel('Voltage (V)')
ax1.set_ylabel(r'Heat ($\frac{k_B T}{cycle}$)')
ax1.plot(voltage,Qdot_x, color='red', marker='.', linewidth=1, markersize=3, label=r'$\dot{Q}_X$')
ax1.plot(voltage,Qdot_y, color='blue', marker='.', linewidth=1, markersize=3, label=r'$\dot{Q}_Y$')
ax1.legend(loc="upper right")
#plt.savefig('results/voltage_sub_heat.png', format='png', bbox_inches='tight', dpi=300)

#fig, ax = plt.subplots(1,1, layout='tight')
ax2.set_title("Voltage v. Subsystems Work")
ax2.grid()
ax2.set_xlabel('Voltage (V)')
ax2.set_ylabel(r'Work ($\frac{k_B T}{cycle}$)')
ax2.plot(voltage,Wdot_x, color='red', marker='.', linewidth=1, markersize=3, label=r'$\dot{W}_X$')
ax2.plot(voltage,Wdot_y, color='blue', marker='.', linewidth=1, markersize=3, label=r'$\dot{W}_Y$')
ax2.legend(loc="upper right")
plt.savefig('results/voltage_sub.png', format='png', bbox_inches='tight', dpi=300)


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
