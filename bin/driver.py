#!/usr/bin/python3

# Program to read transition rates and pass them to main program
from pandas import read_csv
import subprocess

main_cpp = "src/19_states.cpp"
main_exe = "bin/19_states.x"
flags_cpp = ['--std=c++20', '-Wall']

transition_rates = read_csv("data/transition_rates.csv", skiprows=0, sep=',')
parameters = read_csv("data/parameters.csv", skiprows=0, sep=',')
params=[str(i) for i in transition_rates.iloc[:,1]]+[str(i) for i in parameters.iloc[:,1]]

print('Deleting '+main_exe)
subprocess.run(['rm',main_exe],stderr=subprocess.DEVNULL)

print('Compiling '+main_cpp+'...')
subprocess.run(['g++',*flags_cpp,main_cpp,'-o',main_exe])

print('Running '+main_exe+'...')
try:
    subprocess.run([main_exe,*params])
    print('Results saved in results directory.')
except Exception as e:
    print(e)
    print('Could not run '+main_exe+'. Check runtime errors.')
