import sys
sys.path.append('./src')

import PythonBiosystemFramework as pbf
import matplotlib.pyplot as plt

# Create a BioSystem to simulate.
sys = pbf.StochasticBioSystem()

# Add the constant 'k' with a value 0.05.
# This will be the speed of the chemical reaction.
sys.addConstant('C1', 0.7)
sys.addConstant('C11', 0.5)
sys.addConstant('C2', 2)
sys.addConstant('C21', 0.25)

# Let's define substances: A, B, E
# with initial contration values: 10, 0, 1
dAdt = sys.addCompositor('A', 0)
dBdt = sys.addCompositor('B', 0)
dEdt = sys.addCompositor('E', 0)

# Define chemical reaction 'A + E -k> B + E' rates
# Set initial substance concentrations.
# Substance A decreases by the law '-k * A * E' where A, E are current
# concentration of substances.
# Substance B increases by the law 'k * A * E' (as much as A decreases).
# Substance E is a catalyst, so, its concentration doesn't change.
# reaction  = pbf.OdePart(
# 'A + E -k> B + E',
# [dAdt, dBdt, dEdt],
# [pbf.Rate('-k * A * E'), pbf.Rate('k * A * E'), pbf.Rate('0')])

reaction2 = pbf.StochasticPart(
    '* -c> 2A',
    [dAdt, dBdt, dEdt],
    [2, 0, 0],
    '1 * C1'
)

reaction21 = pbf.StochasticPart(
    'A -c> *',
    [dAdt, dBdt, dEdt],
    [-1, 0, 0],
    'A * C11'
)

reaction3 = pbf.StochasticPart(
    '* -A-> B',
    [dAdt, dBdt, dEdt],
    [0, 1, 0],
    'A * C2'
)

reaction31 = pbf.StochasticPart(
    'B -c> *',
    [dAdt, dBdt, dEdt],
    [0, -1, 0],
    'B * C21'
)

# Add the reaction to the simulation.
sys.addPart(reaction2)
sys.addPart(reaction3)
sys.addPart(reaction21)
sys.addPart(reaction31)

# Initialise time points and substance concentration values.
T = None
Y = None

pulses = []

# initial condition
pulses.append(pbf.Pulse(0, 'A', 0))

# spike in some A
pulses.append(pbf.Pulse(4, 'A', 0))

pulses.append(pbf.Pulse(6, 'A', 0))

pulses.append(pbf.Pulse(8, 'A', 0))
# stop the simulation at time 500 with this empty string as the state
# variable parameter
pulses.append(pbf.Pulse(10, '', 0))

# Run pulsed simulation
(T, Y) = sys.run_pulses(pulses)

# Plot the simulation data
plt.figure()
plt.margins(0.05, 0.1)
plt.grid(True)
plt.plot(T, Y[:, sys.compositorIndex('A')], label="A")
plt.plot(T, Y[:, sys.compositorIndex('B')], label="B")
#plt.plot(T, Y[:, sys.compositorIndex('E')], label="E")
plt.legend()
plt.xlabel('Time')
plt.ylabel('Molecules')
plt.show()