import sys
sys.path.append('./src')

import PythonBiosystemFramework as pbf
import matplotlib.pyplot as plt
from scipy.stats import *
import numpy as np

# Create a BioSystem to simulate.
sys = pbf.BioSystem()

# Add the constant 'k' with a value 0.05.
# This will be the speed of the chemical reaction.
sys.addConstant('MM', 10) # IFTs count in flagella
sys.addConstant('aa', 10.0e-9) # length of element
sys.addConstant('vv', hmean([2.0e-6, 3.5e-6])) # average speed of movement of
# IFTs parts in flagella.
sys.addConstant('VV', 10e-9) # decay speed of flagella tip

# Let's define substances: A, B, E
# with initial contration values: 10, 0, 1
dNdt = sys.addCompositor('NN', 1900)# 100, 400, 900

# Flagella length diff. equation.
# Modified to model length by number of parts of flaggela.
# Formula L = a * N was used to make substitution in original diff. eq.
#   L - length
#   a - length of the part
#   N - number of parts
reaction = pbf.OdePart(
    'Flagella Length',
    [dNdt],
    [pbf.Rate('(MM * vv / (2 * aa * NN)) - (VV / aa)')])

# Add the reaction to the simulation.
sys.addPart(reaction)

# Initialise time points and substance concentration values.
T = None
Y = None

# Simulate system with provided reactions for 15000 seconds.
(T, Y) = sys.run([0, 15000], 15000)

# save the result
TY = np.concatenate((np.vstack(T), Y), axis=1)
print(TY)
np.savetxt("flagella_N_1900_ode.csv", TY, delimiter=';')

# T - time points of the simulation.
# Y - a matrix, rows shows the substance concentrations at particular time
# point, columns - substance concentrations change in time.

# Plot the simulation data.
plt.figure()
# Create a 5% (0.05) and 10% (0.1) padding in the
# x and y directions respectively.
plt.margins(0.05, 0.1)
plt.grid(True)
plt.plot(T, Y[:, sys.compositorIndex('NN')], label="NN")
plt.legend()
plt.xlabel('Time')
plt.ylabel('Flagella length in parts')
plt.show()