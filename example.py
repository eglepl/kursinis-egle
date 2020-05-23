import sys
sys.path.append('src')

import PythonBiosystemFramework as pbf
import matplotlib.pyplot as plt

print(pbf)

# Create a BioSystem to simulate.
sys = pbf.BioSystem()

# Add the constant 'k' with a value 0.05.
# This will be the speed of the chemical reaction.
sys.addConstant('k', 0.05)

# Let's define substances: A, B, E
# with initial contration values: 10, 0, 1
dAdt = sys.addCompositor('A', 10)
dBdt = sys.addCompositor('B', 0)
dEdt = sys.addCompositor('E', 1)

# Define chemical reaction 'A + E -k> B + E' rates
# Set initial substance concentrations.
# Substance A decreases by the law '-k * A * E' where A, E are current
# concentration of substances.
# Substance B increases by the law 'k * A * E' (as much as A decreases).
# Substance E is a catalyst, so, its concentration doesn't change.
reaction  = pbf.OdePart(
'A + E -k> B + E',
[dAdt, dBdt, dEdt],
[pbf.Rate('-k * A * E'), pbf.Rate('k * A * E'), pbf.Rate('0')])

# reaction2 = StochasticPart(example.py
#     'A + A -c> 2B',
#     [dAdt, dBdt, dEdt],
#     [-2, 2, 0],
#     'A * (A - 1) / 2.0'
#     'C1'
# )

# Add the reaction to the simulation.
sys.addPart(reaction)

# Initialise time points and substance concentration values.
T = None
Y = None

# Simulate system with provided reactions for 25 seconds.
(T, Y) = sys.run([0, 25], 1000)

# T - time points of the simulation.
# Y - a matrix, rows shows the substance concentrations at particular time
# point, columns - substance concentrations change in time.

# Plot the simulation data.
plt.figure()
plt.plot(T, Y[:, sys.compositorIndex('A')], label="A")
plt.plot(T, Y[:, sys.compositorIndex('B')], label="B")
plt.plot(T, Y[:, sys.compositorIndex('E')], label="E")
plt.legend()
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.show()