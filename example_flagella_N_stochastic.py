import sys
sys.path.append('./src')

import PythonBiosystemFramework.stochastic as pbfs
import matplotlib.pyplot as plt
from scipy.stats import *
import numpy as np

# Create a BioSystem to simulate.
sys = pbfs.BioSystem()

#####################################################################
## Define constants in the system
sys.addConstant('a', 10e-9) # Flagellar unit length, m
sys.addConstant('M', 10) # Number of IFT particles in model
sys.addConstant('lmb_p', 200) # Mean speed of anterograde motion of IFT
# particle, units/s
sys.addConstant('lmb_m', 350) # Mean speed of retrograde motion of IFT
# particle, units/s
sys.addConstant('mu', 1) # mean speed of disassembly in units, s^-1
sys.addConstant('Lavg', 12.727e-6) # Average length of Flagellar, m
sys.addConstant('Navg', 1272.7) # Average number of units in Flagellar, #

sys.addConstant('IFT_init_pos', 0) # Initial value for all IFT state
sys.addConstant('NN', 1900) # Initial Flagellar length in units

#####################################################################
## Define initial states

# IFT position compositors vector
x = {}
for i in range(0, sys.getConstantValue('M')):
    c_name = 'ift' + str(i)
    c = sys.addCompositor(c_name, sys.getConstantValue('IFT_init_pos'))
    x[c_name] = c

c_n = sys.addCompositor('n', sys.getConstantValue('NN'))


#####################################################################
## Define events of the system and its rates

# Degradation event
class DegradationPart(pbfs.Part):

    def __init__(self, name, compositor_name):
        pbfs.Part.__init__(self, name)
        self.compositor_name = compositor_name
        self.mu = None

    def prepare(self, bio_system):
        self.mu = bio_system.getConstantValue('mu')

    def process(self, y, t):
        yy = list(y)
        for j in range(0, len(y)-1):
            if yy[j] == yy[-1]:
                yy[j] -= 1
            elif yy[j] == -yy[-1]:
                yy[j] += 1
        yy[-1] -= 1
        return yy

    def get_rate(self, y, t):
        return self.mu

# IFT move event
class IFTMovePart(pbfs.Part):

    def __init__(self, name, compositor_name):
        pbfs.Part.__init__(self, name)
        self.compositor_name = compositor_name
        self.compositor_index = None
        self.lmb_p = None
        self.lmb_m = None

    def prepare(self, bio_system):
        self.lmb_p = bio_system.getConstantValue('lmb_p')
        self.lmb_m = bio_system.getConstantValue('lmb_m')
        self.compositor_index = bio_system.compositorIndex(self.compositor_name)

    def process(self, y, t):
        yy = list(y)
        xi = yy[self.compositor_index]
        if xi == yy[-1]:
            yy[-1] += 1
            yy[self.compositor_index] = -yy[-1]
        else:
            yy[self.compositor_index] += 1

        return yy

    def get_rate(self, y, t):
        xi = y[self.compositor_index]
        if (0 <= xi) and (xi <= y[-1]):
            return self.lmb_p
        else:
            return self.lmb_m


IFT_move_parts = {}
for i in range(0, sys.getConstantValue('M')):
    ift_compositor_name = "ift" + str(i)
    ift_part_name = "ift" + str(i) + "_move_part"
    ift_part = IFTMovePart(ift_part_name, ift_compositor_name)
    sys.addPart(ift_part)



# Add dissasembly-degradation event
disassembly_part = DegradationPart('degradation_part', 'n')
sys.addPart(disassembly_part)

# Initialise time points and substance concentration values.
T = None
Y = None

# Simulate system with provided reactions for 15000 seconds.
(T, Y) = sys.run([0, 5000])

# save the result
TY = np.concatenate((np.vstack(T), Y), axis=1)
#print(TY)
np.savetxt("flagella_N_1900_stoch_v1.csv", TY, delimiter=';')

# T - time points of the simulation.
# Y - a matrix, rows shows the substance concentrations at particular time
# point, columns - substance concentrations change in time.

# Plot the simulation data.
plt.figure()
# Create a 5% (0.05) and 10% (0.1) padding in the
# x and y directions respectively.
plt.margins(0.05, 0.1)
plt.grid(True)
plt.plot(T, Y[:, sys.compositorIndex('n')], label="Flagellar length in unints, N")
plt.legend()
plt.xlabel('Time')
plt.ylabel('Flagella length in parts')
plt.show()