# -*- coding: utf-8 -*-

import sys
sys.path.append('./src')
sys.path.append('.')

from CellFlagellaOdeBiosystemMu2 import * 
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
    sys.addConstant('_N', 10) # IFT count
    sys.addConstant('_beta', 0.5) # faction of A zone elements that one IFT takes

    sys.addConstant('_delta', 1) # tip length
    sys.addConstant('_eta', 1) # flagella cross-section area
    sys.addConstant('_gamma', 1) # amount of kinesin that one IFT carries
    sys.addConstant('_kba', 0.1) # kinesin transfer rate Zone B to Zone A

    sys.addConstant('_lmb_1', 2) # Mean speed of anterograde motion of IFT
    # particle, units/s
    sys.addConstant('_lmb_p', 2) # Mean speed of anterograde motion of IFT
    # particle, units/s
    sys.addConstant('_lmb_2', 2) # Mean speed of anterograde motion of IFT
    # particle, units/s
    sys.addConstant('_lmb_3', 3.5) # Mean speed of anterograde motion of IFT
    # particle, units/s
    sys.addConstant('_lmb_m', 3.5) # Mean speed of retrograde motion of IFT
    # particle, units/s
    sys.addConstant('_lmb_4', 3.5) # Mean speed of anterograde motion of IFT

    # harmonic mean of all lambda
    sys.addConstant('_omega', 6.0 / (
        1.0/sys.getConstantValue('_lmb_1') + 
        1.0/sys.getConstantValue('_lmb_p') + 
        1.0/sys.getConstantValue('_lmb_2') + 
        1.0/sys.getConstantValue('_lmb_3') + 
        1.0/sys.getConstantValue('_lmb_m') + 
        1.0/sys.getConstantValue('_lmb_4'))
        )

    # particle, units/s
    sys.addConstant('_mu', 0.1) # mean speed of disassembly in units, s^-1


sys = BioSystem_OdeCellFlagella(init_constants)

sys.add_flagella(0, 0, 0)
#sys.add_flagella(1, 0, 0)

sys.setup()

#####################################################################
## Define constants in the system

t_end = 2000

# Initialise time points and substance concentration values.
T = None
Y = None


pulses = [
   pbf.Pulse(0, "Ae", 100),
   pbf.Pulse(0, "Ak", 500),
   pbf.Pulse(300, "Ae", 200),
   pbf.Pulse(800, "Ae", 1),
   pbf.Pulse(t_end, '', 0)
   ]
(T, Y) = sys.run_pulses(pulses)

# Simulate system with provided reactions for 15000 seconds.
#(T, Y) = sys.run([0, t_end])

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

plt.plot(T, Y[:, sys.compositorIndex('Ae')], label="Zone A element consentration")
plt.plot(T, Y[:, sys.compositorIndex('Ak')], label="Zone A kinesin concentration")
for i in range(len(sys.L)):
    plt.plot(T, Y[:, sys.compositorIndex('L' + str(i))], label="Flagellum #" + str(i+1) + " length")
    plt.plot(T, Y[:, sys.compositorIndex('Bk' + str(i))], label="Flagellum #" + str(i+1) + " B kinesin")

plt.legend()
plt.xlabel('Time')
plt.ylabel('Elements')
plt.show()