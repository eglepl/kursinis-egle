# -*- coding: utf-8 -*-

import sys
sys.path.append('./src')
sys.path.append('.')

from CellFlagellaOdeBiosystem import * 
import PythonBiosystemFramework as pbf
import PythonBiosystemFramework.ode as pbfo
import matplotlib.pyplot as plt
from scipy.stats import *
import numpy as np
from datetime import datetime
from random import *
from copy import deepcopy



def print_log(s):
    #print(s)
    None

def init_constants(sys):
    sys.addConstant('_k', 1) 
    u = 1
    sys.addConstant('_alpha', 0.01) # how much IFT attach rate drops on Flagella length increse
    sys.addConstant('_beta', 0.3) # ift container, take fraction
    sys.addConstant('_betab', 0.6 / (u)) # ift container, take fraction

    sys.addConstant('_lmb_1', u*2) # Mean speed of anterograde motion of IFT
    # particle, units/s
    sys.addConstant('_lmb_p', u*2) # Mean speed of anterograde motion of IFT
    # particle, units/s
    sys.addConstant('_lmb_2', u*2) # Mean speed of anterograde motion of IFT
    # particle, units/s
    sys.addConstant('_lmb_3', u*3.5) # Mean speed of anterograde motion of IFT
    # particle, units/s
    sys.addConstant('_lmb_m', u*3.5) # Mean speed of retrograde motion of IFT
    # particle, units/s
    sys.addConstant('_lmb_4', u*3.5) # Mean speed of anterograde motion of IFT
    # particle, units/s
    sys.addConstant('_mu', 1) # mean speed of disassembly in units, s^-1

    sys.addConstant('A_zone_init', 1000) # Initial value for all Flagella state
    # Flagella length per 1 motor (IFT): 1.25 nm
    # decay rate of flagella 0.01 mm/s


sys = BioSystem_OdeCellFlagella(init_constants)

sys.add_flagella(0, 0, 0)

sys.setup()

#####################################################################
## Define constants in the system

t_end = 1000

# Initialise time points and substance concentration values.
T = None
Y = None


#pulses = [
#    pbf.Pulse(0, "L0", 1000),
#    pbf.Pulse(0, "L1", 0),
#    pbf.Pulse(t_end, '', 0)
#    ]
#(T, Y) = sys.run_pulses(pulses)

# Simulate system with provided reactions for 15000 seconds.
(T, Y) = sys.run([0, t_end])

# save the result
TY = np.concatenate((np.vstack(T), Y), axis=1)
#print_log(TY)
#print_log(TY)
#np.savetxt("flagella_N_1900_stoch_0-10_v1.csv", TY, delimiter=';',
# fmt='%10.5f')

filename = "output/" + __file__ + "(" + str(len(sys.L)) + ")"
filename += "_ode_0-" + str(t_end)
filename += "_" + datetime.now().strftime("%Y%m%dT%H%M%S")
filename += "_" + str(randrange(100000))
filename += ".csv"
np.savetxt(filename, TY, delimiter=';', fmt='%10.5f')

# T - time points of the simulation.
# Y - a matrix, rows shows the substance concentrations at particular time
# point, columns - substance concentrations change in time.
#
# Plot the simulation data.

plt.figure()
# Create a 5% (0.05) and 10% (0.1) padding in the
# x and y directions respectively.
plt.margins(0.05, 0.1)
plt.grid(True)

plt.plot(T, Y[:, sys.compositorIndex('A')], label="Zone A")
for i in range(len(sys.L)):
    plt.plot(T, Y[:, sys.compositorIndex('L' + str(i))], label="Flagellum #" + str(i+1) + " length")
    #plt.plot(T, Y[:, sys.compositorIndex('B' + str(i))], label="Flagellum #" + str(i+1) + " B zone")


plt.legend()
plt.xlabel('Time')
plt.ylabel('Elements')
plt.show()