Python BioSystem Framework
==========================

The aim of the project was to create simple, easy accessible and editable framework for synthetic biology research,  based on freely available libraries and programming languages. 

The tool was inspired by MIT university course „20.305x Principles of Synthetic Biology“ provided tool „Part-compositor framework“ which is based on „MatLab“ framework. 

The too has been updated to support discrete stochastic model simulation using Gillespie algorithm.

Python programming language and its non-standard libraries: Sympy, Numpy, and Scipy were used to implement the goal, due to the similarity of the MatLab functionality required.
The created tool can simulate concentrations of substances in time using chemical reaction differential equitations with the specified initial concentration conditions of substances.

The implementation is available on the public github webpage: 
https://github.com/eglepl/pybiosystem_framwork


----------


Requirements
-------------

- Linux (might work with other OS)
- Python 2.7
- Python libraries:
	- SymPy v1.0
	- NumPy v1.11.1
	- SciPy v0.18.1
	- Matplotlib v1.5.3 (data plotting)

----------

Documentation
-------------

See the docs/html/index.php file for documentation reference and
examples.

The PDF version can be found in docs/latex/refman.pdf 

----------

Usage Ode
---------

```
# Add ./src directory to PYTHONPATH

import PythonBiosystemFramework as pbf
import PythonBiosystemFramework.ode as pbfo

import matplotlib.pyplot as plt

# Create a BioSystem to simulate.
sys = pbfo.BioSystem()

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
reaction  = pbfo.Part(
'A + E -k> B + E',
[dAdt, dBdt, dEdt],
[pbfo.Rate('-k * A * E'), pbfo.Rate('k * A * E'), pbfo.Rate('0')])

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
```

Usage Stochastic
----------------

```
import PythonBiosystemFramework as pbf
import PythonBiosystemFramework.stochastic as pbfs
import matplotlib.pyplot as plt

# Create a BioSystem to simulate.
sys = pbfs.BioSystem()

# Add the constant 'k' with a value 0.05.
# This will be the speed of the chemical reaction.
sys.addConstant('C1', 1)
sys.addConstant('C11', 0.5)
sys.addConstant('C2', 2)
sys.addConstant('C21', 0.1)

# Let's define substances: A, B, E
# with initial contration values: 10, 0, 1
dAdt = sys.addCompositor('A', 0)
dBdt = sys.addCompositor('B', 0)
dEdt = sys.addCompositor('E', 0)

reaction2 = pbfs.SimplePart(
    '* -c> 2A', # description
    [dAdt, dBdt, dEdt], # used compositors
    [2, 0, 0], # state change
    '1 * C1' # rate of this part/reaction
)

reaction21 = pbfs.SimplePart(
    'A -c> *',
    [dAdt, dBdt, dEdt],
    [-1, 0, 0],
    'A * C11'
)

reaction3 = pbfs.SimplePart(
    '* -A-> B',
    [dAdt, dBdt, dEdt],
    [0, 1, 0],
    'A * C2'
)

reaction31 = pbfs.SimplePart(
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

# Simulate system with provided reactions for 25 seconds.
(T, Y) = sys.run([0, 25])

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
#plt.plot(T, Y[:, sys.compositorIndex('E')], label="E")
plt.legend()
plt.xlabel('Time')
plt.ylabel('Molecules')
plt.show()
```