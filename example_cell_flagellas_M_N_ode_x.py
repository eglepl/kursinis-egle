import sys
sys.path.append('./src')
sys.path.append('.')

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
    sys.addConstant('k', 1) 
    sys.addConstant('_alpha', 0.001) # how much IFT attach rate drops on Flagella length increse
    sys.addConstant('_beta', 0.01) # ift container, take fraction

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
    # particle, units/s
    sys.addConstant('_mu', 0.5) # mean speed of disassembly in units, s^-1

    sys.addConstant('_kbc', 0.01)
    sys.addConstant('_kck', 1)

    sys.addConstant('A_zone_init', 200) # Initial value for all Flagella state
    # Flagella length per 1 motor (IFT): 1.25 nm
    # decay rate of flagella 0.01 mm/s


sys = pbfo.BioSystem()

init_constants(sys)

# A zone element count
A = sys.addCompositor('A', sys.getConstantValue('A_zone_init'))

# Flagella initial state
L = list()
B = list()
BK = list()
CK = list()

def add_flagella(idx, size, b, bk=0, ck=0):
    c_name = 'L' + str(idx)
    c = sys.addCompositor(c_name, size) # Initial value for all Flagella state
    L.append(c)
    c_name = 'B' + str(idx)
    c = sys.addCompositor(c_name, b) # Initial value for all Flagella state
    B.append(c)
    c_name = 'BK' + str(idx)
    c = sys.addCompositor(c_name, bk) # Initial value for all Flagella state
    BK.append(c)
    c_name = 'CK' + str(idx)
    c = sys.addCompositor(c_name, ck) # Initial value for all Flagella state
    CK.append(c)

add_flagella(0, 100, 0)
#add_flagella(1, 0, 0)


for j in range(0, len(L)):
    reaction  = pbfo.Part(
    'L increase',
    [L[j], A, BK[j]],
    [pbfo.Rate('_beta  * A * _lmb_p'),
    pbfo.Rate('-_beta * 3 * A * _lmb_p'),
    pbfo.Rate('_beta * 3 * A * _lmb_p')])
    sys.addPart(reaction)

for j in range(0, len(L)):
    reaction  = pbfo.Part(
    'Tip decompose',
    [L[j], B[j]],
    [pbfo.Rate('-_mu * (1 + BK' + str(j) + ')'), pbfo.Rate('_mu * (1 + BK' + str(j) + ')')])
    sys.addPart(reaction)

for j in range(0, len(L)):
    reaction  = pbfo.Part(
    'B decrease, A increase',
    [B[j], A],
    [pbfo.Rate('-1 * _beta * B' + str(j) + ' * _lmb_p'), 
    pbfo.Rate('1 * _beta * B' + str(j) + ' * _lmb_p')])
    sys.addPart(reaction)

for j in range(0, len(L)):
    reaction  = pbfo.Part(
    'BK decrease, CK increase',
    [BK[j], CK[j]],
    [pbfo.Rate('-_kbc * BK' + str(j) + ' / (1 + CK' + str(j) + ' / (L' +  str(j) + ' + 1))'),
    pbfo.Rate('_kbc * BK' + str(j) + ' / (1 + CK' + str(j) + ' / (L' +  str(j) + ' + 1))')])
    sys.addPart(reaction)

for j in range(0, len(L)):
    reaction  = pbfo.Part(
    'CK decrease',
    [CK[j], A],
    [pbfo.Rate('-_kck * CK' + str(j)),
    pbfo.Rate('_kck * CK' + str(j))])
    sys.addPart(reaction)

#####################################################################
## Define constants in the system

t_end = 300

# Initialise time points and substance concentration values.
T = None
Y = None


# pulses = [
#     pbf.Pulse(0, "a_zone", 0),
#     pbf.Pulse(0, "flagella_0", Flagella(1000, 0)),
#     pbf.Pulse(0, "flagella_1", Flagella(1000, 0)),
#     pbf.Pulse(50, "flagella_0", Flagella(1, 0)),
#     pbf.Pulse(80, "a_zone", 300),
#     pbf.Pulse(150, '', 0)
#     ]

# (T, Y) = sys.run_pulses(pulses)

# Simulate system with provided reactions for 15000 seconds.
(T, Y) = sys.run([0, t_end])

f0_size_y = [y[sys.compositorIndex('L0')] for y in Y]
f0_b_y = [y[sys.compositorIndex('B0')] for y in Y]

# f1_size_y = [y[sys.compositorIndex('flagella_1')].size for y in Y]
# f1_b_y = [y[sys.compositorIndex('flagella_1')].b for y in Y]

a_zone_y = [y[sys.compositorIndex('A')] for y in Y]

# save the result
TY = np.concatenate((np.vstack(T), Y), axis=1)
#print_log(TY)
#print_log(TY)
#np.savetxt("flagella_N_1900_stoch_0-10_v1.csv", TY, delimiter=';',
# fmt='%10.5f')

filename = "output/cell_flagellas_M_N_(" + str(len(L)) + ")"
filename += "_ode_0-" + str(t_end)
filename += "_" + datetime.now().strftime("%Y%m%dT%H%M%S")
filename += "_" + str(randrange(100000))
filename += ".csv"
#np.savetxt(filename, TY, delimiter=';', fmt='%10.5f')

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

plt.plot(T, Y[:, sys.compositorIndex('A')], label="A zone")
for i in range(len(L)):
    plt.plot(T, Y[:, sys.compositorIndex('L' + str(i))], label="Žiuželio nr. " + str(i) + " ilgis")
    plt.plot(T, Y[:, sys.compositorIndex('B' + str(i))], label="Žiuželio nr. " + str(i) + " B zona")
    plt.plot(T, Y[:, sys.compositorIndex('BK' + str(i))], label="BK" + str(i) + "")
    plt.plot(T, Y[:, sys.compositorIndex('CK' + str(i))], label="CK" + str(i) + "")

# plt.plot(T, f1_size_y, label="Flagellar2 length in unints, N")
# plt.plot(T, f1_b_y, label="Flagellar2 B zone elements count, N")

#plt.plot(T, a_zone_y, label="A zone elements count")

plt.legend()
plt.xlabel('Time')
plt.ylabel('Flagella length in parts')
plt.show()