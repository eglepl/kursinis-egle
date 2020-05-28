import sys
sys.path.append('./src')
sys.path.append('.')

from CellFlagellaStochasticBiosystemNoB import * 
import matplotlib.pyplot as plt
from scipy.stats import *
import numpy as np
import pandas as pd
from utils import *

def print_log(s):
    #print(s)
    None

def init_constants(sys):
    sys.addConstant('_N', 20) # Number of IFT particles in model
    sys.addConstant('_M', 2) # Flagella count
    sys.addConstant('C_max', 10) # ift container, max elements

    sys.addConstant('_alpha', 0.01) # how much IFT attach rate drops on Flagella length increse
    sys.addConstant('_beta', 0.3) # ift container, take fraction
    sys.addConstant('_gamma', 1) # length per element
    sys.addConstant('_zeta', 1)  #elements to dissociate from flagella end.
    sys.addConstant('_phi', 1) # length of one element at flagellum 

    sys.addConstant('lmb_1', 2) # Mean speed of anterograde motion of IFT
    # particle, units/s
    sys.addConstant('lmb_p', 2) # Mean speed of anterograde motion of IFT
    # particle, units/s
    sys.addConstant('lmb_2', 2) # Mean speed of anterograde motion of IFT
    # particle, units/s
    sys.addConstant('lmb_3', 3.5) # Mean speed of anterograde motion of IFT
    # particle, units/s
    sys.addConstant('lmb_m', 3.5) # Mean speed of retrograde motion of IFT
    # particle, units/s
    sys.addConstant('lmb_4', 3.5) # Mean speed of anterograde motion of IFT
    # particle, units/s
    sys.addConstant('_mu', 1) # mean speed of disassembly in units, s^-1

    sys.addConstant('s_p', 1) # 8nm # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2186253/
    sys.addConstant('s_m', 1) # 8.3nm # https://www.nature.com/articles/s41598-018-34549-7
    # Flagella length per 1 motor (IFT): 1.25 nm
    # decay rate of flagella 0.01 mm/s


# Create a BioSystem to simulate.
sys = BioSystem_StochasticCellFlagella(init_constants)

#####################################################################
## Define constants in the system

#t_end = 1000
t_end = 2000

# Initialise time points and substance concentration values.
T = None
Y = None

pulses = [
    pbf.Pulse(0, "A", 1000),
    pbf.Pulse(0, "flagella_1_U", 600),
    pbf.Pulse(0, "flagella_0_U", 0),
    pbf.Pulse(t_end, '', 0)
     ]

(T, Y) = sys.run_pulses(pulses)

# save the result
TY = np.concatenate((np.vstack(T), Y), axis=1)

save_csv(__file__, sys, TY)

# T - time points of the simulation.
# Y - a matrix, rows shows the substance concentrations at particular time
# point, columns - substance concentrations change in time.
#
# Plot the simulation data.

plt.figure(num=__file__)
# Create a 5% (0.05) and 10% (0.1) padding in the
# x and y directions respectively.
plt.margins(0.05, 0.1)
plt.grid(True)

plt.plot(T, Y[:, sys.compositorIndex('A')], label="A zone")
#plt.plot(T, Y[:, sys.compositorIndex('n')], label="Flagellar length in unints, N")
plt.plot(T, Y[:, sys.compositorIndex('flagella_0_U')], label="Flagellar #1")
#plt.plot(T, Y[:, sys.compositorIndex('flagella_0_B')], label="Flagellar #1 B zone")
plt.plot(T, Y[:, sys.compositorIndex('flagella_1_U')], label="Flagellar #2")
#plt.plot(T, Y[:, sys.compositorIndex('flagella_1_B')], label="Flagellar #2 B zone")

# plt.plot(T, f1_size_y, label="Flagellar2 length in unints, N")
# plt.plot(T, f1_b_y, label="Flagellar2 B zone elements count, N")


plt.legend()
plt.xlabel('Time')
plt.ylabel('Elements')
plt.show()
