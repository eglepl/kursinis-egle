import sys
sys.path.append('./src')

import PythonBiosystemFramework as pbf
import matplotlib.pyplot as plt

sys = pbf.BioSystem()
sys.addConstant('k', 0.05)
dAdt = sys.addCompositor('A', 10)
dBdt = sys.addCompositor('B', 0)
dEdt = sys.addCompositor('E', 1)
reaction  = pbf.OdePart(
'A + E -k> B + E',
[dAdt, dBdt, dEdt],
[pbf.Rate('-k * A * E'), pbf.Rate('k * A * E'), pbf.Rate('0')])
sys.addPart(reaction)
T = None
Y = None

pulses = []

# initial condition
pulses.append(pbf.Pulse(0, 'A', 10))
pulses.append(pbf.Pulse(1, 'B', 5))

# spike in some A
pulses.append(pbf.Pulse(100, 'A', 20))

# spike in a bit less A
pulses.append(pbf.Pulse(150, 'A', 5))

# spike in more A again
pulses.append(pbf.Pulse(250, 'A', 10))

# stop the simulation at time 500 with this empty string as the state
# variable parameter
pulses.append(pbf.Pulse(500, '', 0))

# Run pulsed simulation
(T, Y) = sys.run_pulses(pulses, 100)

# Plot the simulation data
plt.figure()
plt.plot(T, Y[:, sys.compositorIndex('A')], label="A")
plt.plot(T, Y[:, sys.compositorIndex('B')], label="B")
plt.plot(T, Y[:, sys.compositorIndex('E')], label="E")
plt.legend()
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.show()
