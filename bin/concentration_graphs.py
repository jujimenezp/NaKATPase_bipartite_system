#!/usr/bin/python3

from pandas import read_csv
from matplotlib import pyplot as plt

data = read_csv("results/proportion_range.dat", sep="\t")
data_Nai = read_csv("results/proportion_range_Nai.dat", sep="\t")
data_Nai_less = read_csv("results/proportion_range_Nai_less.dat", sep="\t")
data_Nao = read_csv("results/proportion_range_Nao.dat", sep="\t")
data_Nao_less = read_csv("results/proportion_range_Nao_less.dat", sep="\t")
data_Ko = read_csv("results/proportion_range_Ko.dat", sep="\t")
data_Ko_less = read_csv("results/proportion_range_Ko_less.dat", sep="\t")
data_Ki = read_csv("results/proportion_range_Ki.dat", sep="\t")
data_Ki_less = read_csv("results/proportion_range_Ki_less.dat", sep="\t")

prop = data.iloc[:,0]
turnover = data.iloc[:,1]
turnover_Nai = data_Nai.iloc[:,1]
turnover_Nai_less = data_Nai_less.iloc[:,1]
turnover_Nao = data_Nao.iloc[:,1]
turnover_Nao_less = data_Nao_less.iloc[:,1]
turnover_Ko = data_Ko.iloc[:,1]
turnover_Ko_less = data_Ko_less.iloc[:,1]
turnover_Ki = data_Ki.iloc[:,1]
turnover_Ki_less = data_Ki_less.iloc[:,1]

fig, ax = plt.subplots(1,1, figsize=(10,4), layout='tight')
ax.set_xscale('log')
ax.set_xlabel(r"$\frac{W_{0,1}}{W_{7,8}}$")
ax.set_ylabel('Probability current')
ax.plot(prop, turnover, color='black')
ax.fill_between(prop,turnover_Nai,turnover_Nai_less, color='red', alpha=0.5, label=r"(Na$^+$)$_{in}$")
ax.fill_between(prop,turnover_Nao_less,turnover_Nao, color='orange', alpha=0.5, label=r"(Na$^+$)$_{out}$")
ax.fill_between(prop,turnover_Ko_less,turnover_Ko, color='green', alpha=0.5, label=r"(K$^+$)$_{out}$")
ax.fill_between(prop,turnover_Ki_less,turnover_Ki, color='blue', alpha=0.5, label=r'(K$^+$)$_{in}$')
ax.legend(loc='upper right')
plt.savefig('results/concentrations_variation.png', format='png', bbox_inches='tight', dpi=300)
