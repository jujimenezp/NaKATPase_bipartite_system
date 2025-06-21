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
base_current = 0.312442 #Taken from the data when change in concentrations is 0

percentage = data_NaO.iloc[:,0]
current_NaO = data_NaO.iloc[:,2]/base_current-1
current_NaI = data_NaI.iloc[:,2]/base_current-1
current_KO = data_KO.iloc[:,2]/base_current-1
current_KI = data_KI.iloc[:,2]/base_current-1

fig, ax = plt.subplots(1,1, layout = 'tight')
ax.set_box_aspect(0.6)
#ax.grid()
ax.set_xlim(-0.2,0.2)
#ax.set_ylim(-0.5,0.75)
#ax.ticklabel_format(axis='y', style='scientific', scilimits=(-4,-4), useMathText=True)
#ax.set_ylim(1e-4,1e-3)
ax.set_xlabel('Concentration change', fontsize=18)
ax.set_ylabel('Probability current change', fontsize=18)
ax.set_xticks([-0.2, -0.15,- 0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2], ['-0.2', '-0.15', '-0.1', '-0.05', '0', '0.05', '0.1', '0.15', '0.2'], fontsize=16)
ax.set_yticks([-0.4, -0.2, 0, 0.2, 0.4, 0.6], ['-0.4', '-0.2', '0', '0.2', '0.4', '0.6'], fontsize=16)
ax.plot(percentage, current_NaO, color='darkblue', linewidth=3, linestyle='dashed', label=r'Na$^+_\text{out}$')
ax.plot(percentage, current_NaI, color='darkblue', linewidth=3, label=r'Na$^+_\text{in}$')
ax.plot(percentage, current_KO, color='#007700', linewidth=3, linestyle='dashed', label=r'K$^+_\text{out}$')
ax.plot(percentage, current_KI, color='#007700', linewidth=3, label=r'K$^+_\text{in}$')
#ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.f%%'))
#ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.f%%'))
#ax.legend(loc='upper center', frameon=False)
ax.text(0.15, 0.57, r'Na$^+_\text{out}$', transform=ax.transAxes, fontsize=16, verticalalignment='top', backgroundcolor='None', color='darkblue')
ax.text(0.73, 0.54, r'K$^+_\text{out}$', transform=ax.transAxes, fontsize=16, verticalalignment='top', backgroundcolor='None', color='#007700')
ax.text(0.18, 0.78, r'K$^+_\text{in}$', transform=ax.transAxes, rotation=-45, transform_rotates_text=True, fontsize=16, verticalalignment='top', backgroundcolor='None', color='#007700')
ax.text(0.67, 0.78, r'Na$^+_\text{in}$', transform=ax.transAxes, rotation=50, transform_rotates_text=True, fontsize=16, verticalalignment='top', backgroundcolor='None', color='darkblue')
plt.savefig('results/concentration_range.png', format='png', bbox_inches='tight', dpi=450)
