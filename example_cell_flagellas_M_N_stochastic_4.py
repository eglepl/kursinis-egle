import sys
sys.path.append('./src')
sys.path.append('.')

from CellFlagellaStochasticBiosystem import * 
import matplotlib.pyplot as plt
from scipy.stats import *
import numpy as np

def print_log(s):
    #print(s)
    None

def init_constants(sys):
    sys.addConstant('alpha', 0.01) # how much IFT attach rate drops on Flagella length increse
    sys.addConstant('beta', 0.1) # ift container, take fraction
    sys.addConstant('betab', 1) # ift container, take fraction
    sys.addConstant('M_IFT', 1000) # ift container, max elements
    sys.addConstant('_gamma', 1) # length per element
    sys.addConstant('zeta', 1)  #elements to dissociate from flagella end.
    g = 5
    sys.addConstant('_N', 10) # Number of IFT particles in model
    sys.addConstant('_M', 1) # Flagella count
    sys.addConstant('lmb_1', g * 2) # Mean speed of anterograde motion of IFT
    # particle, units/s
    sys.addConstant('lmb_p', g * 2) # Mean speed of anterograde motion of IFT
    # particle, units/s
    sys.addConstant('lmb_2', g * 2) # Mean speed of anterograde motion of IFT
    # particle, units/s
    sys.addConstant('lmb_3', g * 3.5) # Mean speed of anterograde motion of IFT
    # particle, units/s
    sys.addConstant('lmb_m', g * 3.5) # Mean speed of retrograde motion of IFT
    # particle, units/s
    sys.addConstant('lmb_4', g * 3.5) # Mean speed of anterograde motion of IFT
    # particle, units/s
    sys.addConstant('mu', 1) # mean speed of disassembly in units, s^-1

    sys.addConstant('s_p', 1) # 8nm # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2186253/
    sys.addConstant('s_m', 1) # 8.3nm # https://www.nature.com/articles/s41598-018-34549-7
    sys.addConstant('A_zone_init', 1000) # Initial value for all Flagella state
    # Flagella length per 1 motor (IFT): 1.25 nm
    # decay rate of flagella 0.01 mm/s


# Create a BioSystem to simulate.
sys = BioSystem_StochasticCellFlagella(init_constants)

#####################################################################
## Define constants in the system

t_end = 1000

# Initialise time points and substance concentration values.
T = None
Y = None

# Simulate system with provided reactions for 15000 seconds.
(T, Y) = sys.run([0, t_end])

# save the result
TY = np.concatenate((np.vstack(T), Y), axis=1)
#print_log(TY)
#print_log(TY)
#np.savetxt("flagella_N_1900_stoch_0-10_v1.csv", TY, delimiter=';',
# fmt='%10.5f')

filename = "output/" + __file__ + "(" + str(sys.getConstantValue('_M')) + "," + str(sys.getConstantValue('_N')) + ")"
filename += "_stoch_0-" + str(t_end)
filename += "_" + datetime.now().strftime("%Y%m%dT%H%M%S")
filename += "_" + str(randrange(100000))
filename += ".csv"
np.savetxt(filename, TY, delimiter=';', fmt='%10.5f')


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

plt.plot(T, Y[:, sys.compositorIndex('a_zone')], label="A zone")
plt.plot(T, Y[:, sys.compositorIndex('flagella_0_L')], label="Flagellar length")
#plt.plot(T, Y[:, sys.compositorIndex('flagella_0_B')], label="Flagellar B zone")

#plt.plot(T, f1_size_y, label="Flagellar2 length in unints, N")
# plt.plot(T, f1_b_y, label="Flagellar2 B zone elements count, N")

plt.legend()
plt.xlabel('Time')
plt.ylabel('Elements')
plt.show()