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

percentage = data_NaO.iloc[:,0]*100
current_NaO = data_NaO.iloc[:,2]
current_NaI = data_NaI.iloc[:,2]
current_KO = data_KO.iloc[:,2]
current_KI = data_KI.iloc[:,2]

fig, ax = plt.subplots(1,1, layout = 'tight')
ax.set_box_aspect(1)
ax.grid()
ax.set_xlim(-20,20)
ax.set_xlabel('Change in base concentration', fontsize=15)
ax.set_ylabel('Probability current (1/s)', fontsize=15)
ax.plot(percentage, current_NaO, color='orange', linewidth=2, label=r'Na$^+_\text{out}$')
ax.plot(percentage, current_NaI, color='red', linewidth=2, label=r'Na$^+_\text{in}$')
ax.plot(percentage, current_KO, color='green', linewidth=2, label=r'K$^+_\text{out}$')
ax.plot(percentage, current_KI, color='blue', linewidth=2, label=r'K$^+_\text{in}$')
ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f%%'))
ax.legend(loc='upper center')
plt.savefig('results/concentration_range.png', format='png', bbox_inches='tight', dpi=450)
#plt.show()
