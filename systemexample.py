from Biosystem import *
from Part import *
from Rate import *
from Pulse import *
import matplotlib.pyplot as plt

# Create a BioSystem to simulate.
sys = BioSystem()

# Add the constant 'k' with a value 0.05.
# This will be the speed of the chemical reaction.
sys.addConstant('a', 5)
sys.addConstant('b', 0.3)
sys.addConstant('V', 100)


# with initial contration values: 10, 1, 1
dBdt = sys.addCompositor('B', 1) #koncentracija litre


reaction  = Part(
'rates',
[dBdt],
[Rate('B * b * ((V + a * t)^(-1) - t*(V + a * t)^(-2))')])


# Add the reaction to the simulation.
sys.addPart(reaction)

# Initialise time points and substance concentration values.
T = None
Y = None

# Simulate system with provid
# d reactions for 25 seconds.
(T, Y) = sys.run([0, 25])

# T - time points of the simulation.
# Y - a matrix, rows shows the substance concentrations at particular time
# point, columns - substance concentrations change in time.

# Plot the simulation data.
plt.figure()
#plt.plot(T, Y[:, sys.compositorIndex('A')], label="A")
plt.plot(T, Y[:, sys.compositorIndex('B')], label="B")
#plt.plot(T, Y[:, sys.compositorIndex('C')], label="C")


plt.legend()
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.show()