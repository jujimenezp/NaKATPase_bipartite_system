#!/usr/bin/python3

# Script to read transition rates and pass them to main program
from pandas import read_csv
from numpy import linspace
from numpy import logspace
import subprocess
import sys

program = str(sys.argv[1])
main_cpp = program
main_exe = "bin/"+program[4:-4]+".x"
flags_cpp = ['--std=c++20', '-Wall']

transition_rates = read_csv("data/transition_rates.csv", skiprows=0, sep=',')
parameters = read_csv("data/parameters.csv", skiprows=0, sep=',')
params=[str(i) for i in transition_rates.iloc[:,1]]+[str(i) for i in parameters.iloc[:,1]]

# One run of the program
def basic_exe():
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

# Many runs of the program changing transmembrane voltage
def neuron_voltage_range():
    voltage_range = linspace(-87e-3,60e-3,1000)
    print(voltage_range)
    with open('results/voltage_range.dat', 'w') as f:
        print("voltage(V)\tturnover(1/s)\tWork(kBT/cycle)\tHeat(kBT/cycle)\tQdot_x\tQdot_y\tWdot_x\tWdot_y\tIdot_x\tIdot_y\tEfficiency", file=f)

    print('Deleting '+main_exe)
    subprocess.run(['rm',main_exe],stderr=subprocess.DEVNULL)

    print('Compiling '+main_cpp+'...')
    subprocess.run(['g++',*flags_cpp,main_cpp,'-o',main_exe])

    for v in voltage_range:
        params[20] = str(v)
        print(v)
        print('Running '+main_exe+'...')
        try:
            subprocess.run([main_exe,*params])
            print('Results saved in results directory.')
        except Exception as e:
            print(e)
            print('Could not run '+main_exe+'. Check runtime errors.')

#Many runs of the program changing proportion of W_01 to W_78
def prop_range():
    prop_range = logspace(-9,9,300)
    with open('results/proportion_range_Ki_less.dat', 'w') as f:
        print("prop_w01_to_w78\tturnover(1/s)\tWork(kBT/cycle)\tHeat(kBT/cycle)\tQdot_x\tQdot_y\tWdot_x\tWdot_y\tIdot_x\tIdot_y\tEfficiency", file=f)

    print('Deleting '+main_exe)
    subprocess.run(['rm',main_exe],stderr=subprocess.DEVNULL)

    print('Compiling '+main_cpp+'...')
    subprocess.run(['g++',*flags_cpp,main_cpp,'-o',main_exe])

    print('params(31): ',params[31])
    for prop in prop_range:
        params[31] = str(prop)
        print(prop)
        print('Running '+main_exe+'...')
        try:
            subprocess.run([main_exe,*params])
            print('Results saved in results directory.')
        except Exception as e:
            print(e)
            print('Could not run '+main_exe+'. Check runtime errors.')

    print(prop_range)

# Many runs of the program changing the concentration of ions
def concentrations_range():
    base_concentration = float(params[26])
    percentage_change = linspace(-0.2,0.2,101)
    with open('results/concentration_range_KIn.dat', 'w') as f:
        print("percentange_change\tconcentration\tprob_current(1/s)\tWork(kBT/cycle)\tHeat(kBT/cycle)\tQdot_x\tQdot_y\tWdot_x\tWdot_y\tIdot_x\tIdot_y\tEfficiency", file=f)

    print('Deleting '+main_exe)
    subprocess.run(['rm',main_exe],stderr=subprocess.DEVNULL)

    print('Compiling '+main_cpp+'...')
    subprocess.run(['g++',*flags_cpp,main_cpp,'-o',main_exe])

    print('params(26): ', base_concentration)

    with open('results/concentration_range_KIn.dat', 'a') as f:
        for perc in percentage_change:
            params[26] = str(base_concentration*(1+perc))
            print('Running '+main_exe+'...')
            try:
                print(f'{perc:.4f}\t', end='', file=f, flush=True)
                subprocess.run([main_exe,*params])
                print('Results saved in results directory.')
            except Exception as e:
                print(e)
                print('Could not run '+main_exe+'. Check runtime errors.')

if __name__=='__main__':
    neuron_voltage_range()
