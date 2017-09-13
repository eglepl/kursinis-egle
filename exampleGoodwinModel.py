
from Biosystem import *
from Part import *
from Rate import *
from Pulse import *
import matplotlib.pyplot as plt

# Create a BioSystem to simulate.
sys = BioSystem()

sys.addConstant('k1', 0.01)
sys.addConstant('d', 0.1) #constants: k2=k4=k6=d
sys.addConstant('k3', 1)
sys.addConstant('k5', 0.0005)
sys.addConstant('K', 0)

dAdt = sys.addCompositor('A', 1)
dBdt = sys.addCompositor('B', 1)
dCdt = sys.addCompositor('C', 1)
dEdt = sys.addCompositor('E', 1)


reaction  = Part(
'Fist',
[dAdt, dBdt, dCdt, dEdt],
[Rate('k1/(K+C^6)'), Rate('k3 * A'), Rate('k5 * B'), Rate('0')])

# degradation of A to some E
deg_A  = Part(
'A -> E',
[dAdt, dBdt, dCdt, dEdt],
[Rate('-d * A'), Rate('0'), Rate('0'), Rate('d * A')])

# degradation of B to some E
deg_B  = Part(
'B -> E',
[dAdt, dBdt, dCdt, dEdt],
[Rate('0'), Rate('-d * B'), Rate('0'), Rate('d * B')])

# degradation of C to some E
deg_C  = Part(
'C -> E',
[dAdt, dBdt, dCdt, dEdt],
[Rate('0'), Rate('0'), Rate('-d * C'), Rate('d * C')])

# Add the reaction to the simulation.
sys.addPart(reaction)
sys.addPart(deg_A)
sys.addPart(deg_B)
sys.addPart(deg_C)

# Initialise time points and substance concentration values.
T = None
Y = None

# Simulate system with provided reactions for 25 seconds.

(T, Y) = sys.run([0, 1000])

Z = np.concatenate((np.vstack(T), Y), axis=1)
#np.savetxt('goodwin_model_pyhton.csv', Z, fmt='%10.10f')

# T - time points of the simulation.
# Y - a matrix, rows shows the substance concentrations at particular time
# point, columns - substance concentrations change in time.

# Plot the simulation data.
plt.figure(figsize=(10, 4))
plt.grid(True)
plt.plot(T, Y[:, sys.compositorIndex('A')], label="A")
plt.plot(T, Y[:, sys.compositorIndex('B')], label="B")
plt.plot(T, Y[:, sys.compositorIndex('C')], label="C")
plt.legend()
plt.xlabel('Time, s')
plt.ylabel('Concentration, M')
plt.show()
