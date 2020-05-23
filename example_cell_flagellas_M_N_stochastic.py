import sys
sys.path.append('./src')

import PythonBiosystemFramework as pbf
import PythonBiosystemFramework.stochastic as pbfs
import matplotlib.pyplot as plt
from scipy.stats import *
import numpy as np
from datetime import datetime
from random import *
from copy import deepcopy

def print_log(s):
    #print(s)
    None

class IFT:
    def __init__(self, x, d, c, z):
        self.x = x
        self.d = d
        self.c = c
        self.z = z

    def copy(self):
        return IFT(self.x, self.d, self.c, self.z)

    def __str__(self):
        return "{x:" + str(self.x) + ", d:"+ str(self.d) + ", c:" + str(self.c) + ", z:" + str(self.z)  + "}"


class Flagella:
    def __init__(self, size, b):
        self.size = size
        self.b = b

    def copy(self):
        return Flagella(self.size, self.b)

    def __str__(self):
        return "{size:" + str(self.size) + ",  b:"+ str(self.b) + "}"


# Create a BioSystem to simulate.
sys = pbfs.BioSystem()

#####################################################################
## Define constants in the system

t_end = 500

sys.addConstant('alpha', 0.01) # how much IFT attach rate drops on Flagella length increse
sys.addConstant('beta', 0.01) # ift container, take fraction
sys.addConstant('M_IFT', 100) # ift container, max elements
sys.addConstant('_gamma', 1) # length per element
sys.addConstant('zeta', 1)  #elements to dissociate from flagella end.

sys.addConstant('_N', 10) # Number of IFT particles in model
sys.addConstant('_M', 1) # Flagella count
sys.addConstant('lmb_1', 200) # Mean speed of anterograde motion of IFT
# particle, units/s
sys.addConstant('lmb_p', 200) # Mean speed of anterograde motion of IFT
# particle, units/s
sys.addConstant('lmb_2', 200) # Mean speed of anterograde motion of IFT
# particle, units/s
sys.addConstant('lmb_3', 350) # Mean speed of anterograde motion of IFT
# particle, units/s
sys.addConstant('lmb_m', 350) # Mean speed of retrograde motion of IFT
# particle, units/s
sys.addConstant('lmb_4', 350) # Mean speed of anterograde motion of IFT
# particle, units/s
sys.addConstant('mu', 1) # mean speed of disassembly in units, s^-1

sys.addConstant('s_p', 1) # 8nm # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2186253/
sys.addConstant('s_m', 1) # 8.3nm # https://www.nature.com/articles/s41598-018-34549-7
sys.addConstant('A_zone_init', 100000) # Initial value for all Flagella state
# Flagella length per 1 motor (IFT): 1.25 nm
# decay rate of flagella 0.01 mm/s


#####################################################################

## Define initial states

# A zone element count
c_name = 'a_zone'
c = sys.addCompositor(c_name, sys.getConstantValue('A_zone_init'))

# IFT position compositors vector
ifts = {}
for i in range(0, sys.getConstantValue('_N')):
    c_name = 'ift_' + str(i)
    c = sys.addCompositor(c_name, IFT(-1, -1, 0, 'flagella_0')) # Initial value for all IFT state
    ifts[c_name] = c

# Flagella initial state
flagella = {}
for j in range(0, sys.getConstantValue('_M')):
    c_name = 'flagella_' + str(j)
    c = sys.addCompositor(c_name, Flagella(0, 0)) # Initial value for all Flagella state
    flagella[c_name] = c

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
        self.zeta = bio_system.getConstantValue('zeta')
        self.compositor_index = bio_system.compositorIndex(self.compositor_name)
        self.ifts_compositor_indexes = [bio_system.compositorIndex(name) for name in ifts.keys()]

    def process(self, y, t):
        print_log("DegradationPart")
        yy = deepcopy(y)
        flagella = yy[self.compositor_index]
        disociate = min(flagella.size, self.zeta)
        flagella.size -= disociate
        flagella.b += disociate
        if flagella.size <= 0:
            flagella.size = 0
        for ift_index in self.ifts_compositor_indexes:
            if(yy[ift_index].x > flagella.size):
                yy[ift_index].x = flagella.size
        return yy

    def get_rate(self, y, t):
        return self.mu

class IFTZoneAAttachPart(pbfs.Part):
    def __init__(self, part_name, ift_name, flagella_name):
        pbfs.Part.__init__(self, part_name)
        self.ift_name = ift_name
        self.flagella_name = "" + flagella_name

    def prepare(self, bio_system):
        self.lmb_1 = bio_system.getConstantValue('lmb_1')
        self.alpha = bio_system.getConstantValue('alpha')
        self.beta = bio_system.getConstantValue('beta')
        self.m_ift = bio_system.getConstantValue('M_IFT')
        self.ift_index = bio_system.compositorIndex(self.ift_name)
        self.flagella_index = bio_system.compositorIndex(self.flagella_name)
        self.a_zone_index = bio_system.compositorIndex("a_zone")

    def process(self, y, t):
        print_log("IFTZoneAAttachPart: " + self.ift_name + "A->" + self.flagella_name)
        yy = deepcopy(y)
        ift = yy[self.ift_index]
        flagella = yy[self.flagella_index]
        ift.x = 0
        ift.d = +1
        ift.c = max(0, min(self.m_ift, self.beta * yy[self.a_zone_index]))
        yy[self.a_zone_index] -= ift.c
        ift.z = self.flagella_name
        return yy

    def get_rate(self, y, t):
        ift = y[self.ift_index]
        if(ift.x == -1 and ift.d == -1):
            flagella = y[self.flagella_index]
            return float(self.lmb_1) / ( 1 + (self.alpha * flagella.size))
        return 0

class IFTFlagellaPart(pbfs.Part):


    def __init__(self, part_name, ift_name):
        pbfs.Part.__init__(self, part_name)
        self.ift_name = ift_name
        ift_index = None
        bio_system = None

    def prepare(self, bio_system):
        self.ift_index = bio_system.compositorIndex(self.ift_name)
        self.bio_system = bio_system

    def get_ift_flagella(self, y):
        ift = y[self.ift_index]
        flagella_index = self.bio_system.compositorIndex(ift.z)
        return y[flagella_index]


class IFTMovePlus(IFTFlagellaPart):
    def __init__(self, part_name, ift_name):
        IFTFlagellaPart.__init__(self, part_name, ift_name)

    def prepare(self, bio_system):
        IFTFlagellaPart.prepare(self, bio_system)
        self.lmb_p = bio_system.getConstantValue('lmb_p')
        self.s_p = bio_system.getConstantValue('s_p')

    def process(self, y, t):
        print_log("IFTMovePlus" + self.ift_name + "-M->" + y[self.ift_index].z)
        yy = deepcopy(y)
        yy[self.ift_index].x += self.s_p
        print_log(yy[self.ift_index])
        return yy

    def get_rate(self, y, t):
        ift = y[self.ift_index]
        flagella = self.get_ift_flagella(y)
        if(ift.x <= flagella.size - self.s_p and ift.d == +1):
            return float(self.lmb_p)
        return 0

class IFTMoveMinus(IFTFlagellaPart):
    def __init__(self, part_name, ift_name):
        IFTFlagellaPart.__init__(self, part_name, ift_name)

    def prepare(self, bio_system):
        IFTFlagellaPart.prepare(self, bio_system)
        self.lmb_m = bio_system.getConstantValue('lmb_m')
        self.s_m = bio_system.getConstantValue('s_m')

    def process(self, y, t):
        print_log("IFTMoveMinus" + self.ift_name + "<-M-" + y[self.ift_index].z)
        yy = deepcopy(y)
        ift = yy[self.ift_index]
        ift.x -= self.s_m
        print_log(ift)
        return yy

    def get_rate(self, y, t):
        ift = y[self.ift_index]
        flagella = self.get_ift_flagella(y)
        if(ift.x >= self.s_m and ift.d == -1):
                return float(self.lmb_m)
        return 0


class IFTDetachToZoneB(IFTFlagellaPart):
    def __init__(self, part_name, ift_name):
        IFTFlagellaPart.__init__(self, part_name, ift_name)

    def prepare(self, bio_system):
        IFTFlagellaPart.prepare(self, bio_system)
        self.lmb_2 = bio_system.getConstantValue('lmb_2')
        self.s_p = bio_system.getConstantValue('s_p')

    def process(self, y, t):
        print_log("IFTDetachToZoneB" + self.ift_name + "->B" + y[self.ift_index].z)
        yy = deepcopy(y)
        ift = yy[self.ift_index]
        flagella = self.get_ift_flagella(yy)
        ift.x = -1
        ift.d = +1
        flagella.size += ift.c
        ift.c = 0
        return yy

    def get_rate(self, y, t):
        ift = y[self.ift_index]
        flagella = self.get_ift_flagella(y)
        if(ift.x > flagella.size - self.s_p and ift.d == +1):
            return float(self.lmb_2)
        return 0


class IFTZoneBAttachPart(IFTFlagellaPart):
    def __init__(self, part_name, ift_name):
        IFTFlagellaPart.__init__(self, part_name, ift_name)

    def prepare(self, bio_system):
        IFTFlagellaPart.prepare(self, bio_system)
        self.lmb_3 = bio_system.getConstantValue('lmb_3')
        self.beta = bio_system.getConstantValue('beta')
        self.m_ift = bio_system.getConstantValue('M_IFT')

    def process(self, y, t):
        print_log("IFTZoneBAttachPart" + self.ift_name + "<-B" + y[self.ift_index].z)
        yy = deepcopy(y)
        ift = yy[self.ift_index]
        flagella = self.get_ift_flagella(yy)
        ift.x = flagella.size
        ift.d = -1
        elems = max(0, min(self.m_ift, self.beta * flagella.b))
        ift.c = elems
        flagella.b -= elems
        return yy

    def get_rate(self, y, t):
        ift = y[self.ift_index]
        if(ift.x == -1 and ift.d == +1):
            return float(self.lmb_3)
        return 0

class IFTDetachToZoneA(IFTFlagellaPart):
    def __init__(self, part_name, ift_name):
        IFTFlagellaPart.__init__(self, part_name, ift_name)

    def prepare(self, bio_system):
        IFTFlagellaPart.prepare(self, bio_system)
        self.lmb_4 = bio_system.getConstantValue('lmb_4')
        self.s_m = bio_system.getConstantValue("s_m")
        self.a_zone_index = bio_system.compositorIndex("a_zone")

    def process(self, y, t):
        print_log("IFTDetachToZoneA"+ self.ift_name + "A<-" + y[self.ift_index].z)
        yy = deepcopy(y)
        ift = yy[self.ift_index]
        ift.x = -1
        ift.d = -1
        yy[self.a_zone_index] += ift.c
        ift.c = 0
        return yy

    def get_rate(self, y, t):
        ift = y[self.ift_index]
        print_log("XXX")
        print_log(self.s_m)
        print_log(ift.x <= self.s_m)
        if(ift.x <= self.s_m and ift.d == -1):
            return float(self.lmb_4)
        return 0


IFT_move_parts = {}
for i in range(0, sys.getConstantValue('_N')):
    ift_compositor_name = "ift_" + str(i)
    # IFT attach to flagella
    for j in range(0, sys.getConstantValue('_M')):
        flagella_compositor_name = "flagella_" + str(j)
        part_name = ift_compositor_name + "_attach_" + flagella_compositor_name 
        attach_part = IFTZoneAAttachPart(part_name, ift_compositor_name, flagella_compositor_name)
        sys.addPart(attach_part)

    # IFT move plus
    part_name = ift_compositor_name + "_move_plus" 
    move_plus_part = IFTMovePlus(part_name, ift_compositor_name)
    sys.addPart(move_plus_part)

    #IFT detach to zone B
    part_name = ift_compositor_name + "_detach_to_b" 
    detach_to_b_part = IFTDetachToZoneB(part_name, ift_compositor_name)
    sys.addPart(detach_to_b_part)

    #IFT attach to flagella at zone B
    part_name = ift_compositor_name + "_attach_b" 
    attach_zone_b_part = IFTZoneBAttachPart(part_name, ift_compositor_name)
    sys.addPart(attach_zone_b_part)

    # IFT move plus
    part_name = ift_compositor_name + "_move_minus" 
    move_minus_part = IFTMoveMinus(part_name, ift_compositor_name)
    sys.addPart(move_minus_part)

    #IFT detach to zone A
    part_name = ift_compositor_name + "_detach_to_a" 
    detach_to_a_part = IFTDetachToZoneA(part_name, ift_compositor_name)
    sys.addPart(detach_to_a_part)

for i in range(0, sys.getConstantValue('_M')):
    #IFT detach to zone A
    flagella_compositor_name = "flagella_" + str(j)
    part_name = flagella_compositor_name + "_degradation" 
    flagella_degradation_part = DegradationPart(part_name, flagella_compositor_name)
    sys.addPart(flagella_degradation_part)

# Initialise time points and substance concentration values.
T = None
Y = None


# pulses = [
#     pbf.Pulse(0, "a_zone", 0),
#     pbf.Pulse(0, "flagella_0", Flagella(1000, 0)),
#     pbf.Pulse(0, "flagella_1", Flagella(1000, 0)),
#     pbf.Pulse(50, "flagella_0", Flagella(1, 0)),
#     pbf.Pulse(80, "a_zone", 300),
#     pbf.Pulse(150, '', 0)
#     ]

# (T, Y) = sys.run_pulses(pulses)

# Simulate system with provided reactions for 15000 seconds.
(T, Y) = sys.run([0, t_end])

f0_size_y = [y[sys.compositorIndex('flagella_0')].size for y in Y]
f0_b_y = [y[sys.compositorIndex('flagella_0')].b for y in Y]

# f1_size_y = [y[sys.compositorIndex('flagella_1')].size for y in Y]
# f1_b_y = [y[sys.compositorIndex('flagella_1')].b for y in Y]

a_zone_y = [y[sys.compositorIndex('a_zone')] for y in Y]

# save the result
TY = np.concatenate((np.vstack(T), Y), axis=1)
#print_log(TY)
#print_log(TY)
#np.savetxt("flagella_N_1900_stoch_0-10_v1.csv", TY, delimiter=';',
# fmt='%10.5f')

filename = "output/cell_flagellas_M_N_(" + str(sys.getConstantValue('_M')) + "," + str(sys.getConstantValue('_N')) + ")"
filename += "_stoch_0-" + str(t_end)
filename += "_" + datetime.now().strftime("%Y%m%dT%H%M%S")
filename += "_" + str(randrange(100000))
filename += ".csv"
#np.savetxt(filename, TY, delimiter=';', fmt='%10.5f')

# T - time points of the simulation.
# Y - a matrix, rows shows the substance concentrations at particular time
# point, columns - substance concentrations change in time.
#
# Plot the simulation data.

plt.figure()
# Create a 5% (0.05) and 10% (0.1) padding in the
# x and y directions respectively.
plt.margins(0.05, 0.1)
plt.grid(True)

#plt.plot(T, Y[:, sys.compositorIndex('n')], label="Flagellar length in unints, N")
plt.plot(T, f0_size_y, label="Flagellar1 length in unints, N")
plt.plot(T, f0_b_y, label="Flagellar1 B zone elements count, N")

# plt.plot(T, f1_size_y, label="Flagellar2 length in unints, N")
# plt.plot(T, f1_b_y, label="Flagellar2 B zone elements count, N")

#plt.plot(T, a_zone_y, label="A zone elements count")

plt.legend()
plt.xlabel('Time')
plt.ylabel('Flagella length in parts')
plt.show()