
from Biosystem import *
from Part import *
from Rate import *
from Pulse import *
import matplotlib.pyplot as plt

# Create a BioSystem to simulate.
sys = BioSystem()

# Add the constant 'k' with a value 0.05.
# This will be the speed of the chemical reaction.
sys.addConstant('k', 0.1)

# Let's define substances: A, B, E
# with initial contration values: 10, 0, 1
dAdt = sys.addCompositor('A', 10)
dBdt = sys.addCompositor('B', 0)
dRdt = sys.addCompositor('R', 1)

# Define chemical reaction 'A + E -k> B + E' rates
# Set initial substance concentrations.
# Substance A decreases by the law '-k * A * E' where A, E are current
# concentration of substances.
# Substance B increases by the law 'k * A * E' (as much as A decreases).
# Substance E is a catalyst, so, its concentration doesn't change.
reaction  = Part(
'A + R -k> B + R',
[dAdt, dBdt, dRdt],
[Rate('-k * A * R'), Rate('k * A * R'), Rate('0')])

# Add the reaction to the simulation.
sys.addPart(reaction)


# Initialise time points and substance concentration values.
T = None
Y = None

# Simulate system with provided reactions for 25 seconds.
(T, Y) = sys.run([0, 50])


Z = np.concatenate((np.vstack(T), Y), axis=1)
#np.savetxt('simple_enzymatic_reaction_pyhton.csv', Z, fmt='%10.10f')

# T - time points of the simulation.
# Y - a matrix, rows shows the substance concentrations at particular time
# point, columns - substance concentrations change in time.

# Plot the simulation data.
plt.figure(figsize=(10, 4))

plt.plot(T, Y[:, sys.compositorIndex('A')], label="A")
plt.plot(T, Y[:, sys.compositorIndex('B')], label="B")
plt.plot(T, Y[:, sys.compositorIndex('R')], label="R")

plt.grid(True)
plt.legend()
plt.xlabel('Time, s')
plt.ylabel('Concentration, M')
plt.tight_layout()
plt.show()
