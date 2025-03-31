#!/usr/bin/python3

from pandas import read_csv
from matplotlib import pyplot as plt
import matplotlib as mpl
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'

data_NaO = read_csv('results/concentration_range_NaOut.dat', sep="\t")
data_NaI = read_csv('results/concentration_range_NaIn.dat', sep="\t")
data_KO = read_csv('results/concentration_range_KOut.dat', sep="\t")
data_KI = read_csv('results/concentration_range_KIn.dat', sep="\t")
base_current = 0.000430288 #Taken from the data when change in concentrations is 0

percentage = data_NaO.iloc[:,0]
current_NaO = data_NaO.iloc[:,2]/base_current-1
current_NaI = data_NaI.iloc[:,2]/base_current-1
current_KO = data_KO.iloc[:,2]/base_current-1
current_KI = data_KI.iloc[:,2]/base_current-1

fig, ax = plt.subplots(1,1, layout = 'tight')
ax.set_box_aspect(0.6)
#ax.grid()
ax.set_xlim(-0.2,0.2)
ax.set_ylim(-0.5,0.75)
#ax.ticklabel_format(axis='y', style='scientific', scilimits=(-4,-4), useMathText=True)
#ax.set_ylim(1e-4,1e-3)
ax.set_xlabel('Concentration variation', fontsize=18)
ax.set_ylabel('Probability current variation', fontsize=18)
ax.set_xticks([-0.2, -0.15,- 0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2], ['-0.2', '-0.15', '-0.1', '-0.05', '0', '0.05', '0.1', '0.15', '0.2'], fontsize=16)
ax.set_yticks([-0.4, -0.2, 0, 0.2, 0.4, 0.6], ['-0.4', '-0.2', '0', '0.2', '0.4', '0.6'], fontsize=16)
ax.plot(percentage, current_NaO, color='orange', linewidth=5, label=r'Na$^+_\text{out}$')
ax.plot(percentage, current_NaI, color='red', linewidth=2, label=r'Na$^+_\text{in}$')
ax.plot(percentage, current_KO, color='green', linestyle='dashed', linewidth=2, label=r'K$^+_\text{out}$')
ax.plot(percentage, current_KI, color='blue', linewidth=2, label=r'K$^+_\text{in}$')
#ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.f%%'))
#ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.f%%'))
#ax.legend(loc='upper center', frameon=False)
ax.text(0.15, 0.49, r'Na$^+_\text{out}$', transform=ax.transAxes, fontsize=16, verticalalignment='top', backgroundcolor='None', color='orange')
ax.text(0.73, 0.49, r'K$^+_\text{out}$', transform=ax.transAxes, fontsize=16, verticalalignment='top', backgroundcolor='None', color='green')
ax.text(0.18, 0.73, r'K$^+_\text{in}$', transform=ax.transAxes, rotation=-45, transform_rotates_text=True, fontsize=16, verticalalignment='top', backgroundcolor='None', color='blue')
ax.text(0.67, 0.75, r'Na$^+_\text{in}$', transform=ax.transAxes, rotation=53, transform_rotates_text=True, fontsize=16, verticalalignment='top', backgroundcolor='None', color='red')
plt.savefig('results/concentration_range.png', format='png', bbox_inches='tight', dpi=450)
#plt.show()
