import sys
sys.path.append('./src')

import PythonBiosystemFramework as pbf
import matplotlib.pyplot as plt
from scipy.stats import *
import numpy as np

# Create a BioSystem to simulate.
sys = pbf.BioSystem()

aa = 10.0e-9
NN = 100

# Add the constant 'k' with a value 0.05.
# This will be the speed of the chemical reaction.
sys.addConstant('MM', 10) # IFTs count in flagella
sys.addConstant('aa', aa) # length of element
sys.addConstant('vv', hmean([2.0e-6, 3.5e-6])) # average speed of movement of
# IFTs parts in flagella.
sys.addConstant('VV', 10e-9) # decay speed of flagella tip

# Let's define substances: A, B, E
# with initial contration values: 10, 0, 1
dLdt = sys.addCompositor('L', NN * aa)

# Define chemical reaction 'A + E -k> B + E' rates
# Set initial substance concentrations.
# Substance A decreases by the law '-k * A * E' where A, E are current
# concentration of substances.
# Substance B increases by the law 'k * A * E' (as much as A decreases).
# Substance E is a catalyst, so, its concentration doesn't change.
reaction = pbf.OdePart(
'Flagella Length',
[dLdt],
[pbf.Rate('(MM * aa * vv / (2 * L)) - VV')])

# Add the reaction to the simulation.
sys.addPart(reaction)

# Initialise time points and substance concentration values.
T = None
Y = None

# Simulate system with provided reactions for 25 seconds.
(T, Y) = sys.run([0, 15000], 15000)

TY = np.concatenate((np.vstack(T), Y), axis=1)
print(TY)

# T - time points of the simulation.
# Y - a matrix, rows shows the substance concentrations at particular time
# point, columns - substance concentrations change in time.

# Plot the simulation data.
plt.figure()
# Create a 5% (0.05) and 10% (0.1) padding in the
# x and y directions respectively.
plt.margins(0.05, 0.1)
plt.grid(True)
plt.plot(T, Y[:, sys.compositorIndex('A')], label="A")
plt.plot(T, Y[:, sys.compositorIndex('B')], label="B")
plt.plot(T, Y[:, sys.compositorIndex('E')], label="E")
plt.legend()
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.show()