#!/usr/bin/python3

# Program to read transition rates and pass them to main program
from pandas import read_csv
import subprocess

T=300 #Temperature
transition_rates = read_csv("data/transition_rates.csv", header=None, sep='    ', engine='python')
parameters = read_csv("data/parameters.csv", header=None, sep='    ', engine='python')
params=[str(i) for i in transition_rates.iloc[:,1]]+[str(i) for i in parameters.iloc[:,1]]

subprocess.run(['g++','--std=c++20','-Wall','src/main.cpp','-o','bin/main.x'])
subprocess.run(['./bin/main.x',*params])
