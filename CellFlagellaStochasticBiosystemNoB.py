import sys
sys.path.append("./src")

import PythonBiosystemFramework as pbf
import PythonBiosystemFramework.stochastic as pbfs
from scipy.stats import *
import numpy as np
from datetime import datetime
from random import *
from copy import deepcopy
import re

def print_log(s):
    #print(s)
    None

class SubCompositorWrapper:
    def __init__(self, biosystem, compositor_name, index):
        self.biosystem = biosystem
        self.compositor_name = compositor_name
        self.index = index

    def set_y(self, y):
        self.y = y

    def get_subname(self, key):
        name = self.compositor_name + "_" + str(int(self.index)) + "_" + key
        #print("Accessing: '" + name + "'")
        return name

    def __getitem__(self, key):
        name = self.get_subname(key)
        idx = self.biosystem.compositorIndex(name)
        if idx == None:
            print("Name not found: " + name) 
        return self.y[idx]

    def __setitem__(self, key, value):
        idx = self.biosystem.compositorIndex(self.get_subname(key))
        self.y[idx] = value
        return self.y


#####################################################################
## Define events of the system and its rates
# Degradation event
class DegradationPart(pbfs.Part):

    def __init__(self, name, flagella_wrapper):
        pbfs.Part.__init__(self, name)
        self.flagella_wrapper = flagella_wrapper
        self.mu = None

    def prepare(self, bio_system):
        self.mu = bio_system.getConstantValue("_mu")
        self.zeta = bio_system.getConstantValue("_zeta")
        self.phi = bio_system.getConstantValue("_phi")
        self.all_ift_wrappers = bio_system.ift_wrappers
        self.a_zone_index = bio_system.compositorIndex("A")

    def process(self, y, t):
        print_log("DegradationPart")
        yy = list(y)
        self.flagella_wrapper.set_y(yy)

        deconstruct = min(self.flagella_wrapper["U"], self.zeta)
        self.flagella_wrapper["U"] -= deconstruct
        yy[self.a_zone_index] += deconstruct
        for ift_wrapper in self.all_ift_wrappers:
            ift_wrapper.set_y(yy)
            if(ift_wrapper["X"] > self.phi * self.flagella_wrapper["U"]):
                ift_wrapper["X"] = self.phi * self.flagella_wrapper["U"]
        return yy

    def get_rate(self, y, t):
        return self.mu

class IFTZoneAAttachPart(pbfs.Part):
    def __init__(self, part_name, ift_wrapper, flagella_wrapper):
        pbfs.Part.__init__(self, part_name)
        self.ift_wrapper = ift_wrapper
        self.flagella_wrapper = flagella_wrapper

    def prepare(self, bio_system):
        self.lmb_1 = bio_system.getConstantValue("lmb_1")
        self.alpha = bio_system.getConstantValue("_alpha")
        self.beta = bio_system.getConstantValue("_beta")
        self.c_max = bio_system.getConstantValue("C_max")
        self.a_zone_index = bio_system.compositorIndex("A")
        self.phi = bio_system.getConstantValue("_phi")

    def process(self, y, t):
        print_log("IFTZoneAAttachPart: ")
        yy = list(y)
        self.ift_wrapper.set_y(yy)
        self.ift_wrapper["X"] = 0
        self.ift_wrapper["D"] = +1
        self.ift_wrapper["C"] = max(0, min(self.c_max, self.beta * yy[self.a_zone_index]))
        yy[self.a_zone_index] -= self.ift_wrapper["C"]
        self.ift_wrapper["Z"] = self.flagella_wrapper.index
        return yy

    def get_rate(self, y, t):
        self.ift_wrapper.set_y(y)
        self.flagella_wrapper.set_y(y)
        if(self.ift_wrapper["X"] == -1 and self.ift_wrapper["D"] == -1):
            return float(self.lmb_1) / ( 1 + (self.alpha * self.phi * self.flagella_wrapper["U"]))
        return 0

class IFTFlagellaPart(pbfs.Part):

    def __init__(self, part_name, ift_wrapper):
        pbfs.Part.__init__(self, part_name)
        self.ift_wrapper = ift_wrapper

    def prepare(self, bio_system):
        self.bio_system = bio_system


    def get_ift_flagella(self, y):
        self.ift_wrapper.set_y(y)
        flagella_wrapper = SubCompositorWrapper(self.bio_system, "flagella", self.ift_wrapper["Z"])
        flagella_wrapper.set_y(y)
        return flagella_wrapper


class IFTMovePlus(IFTFlagellaPart):
    def __init__(self, part_name, ift_wrapper):
        IFTFlagellaPart.__init__(self, part_name, ift_wrapper)

    def prepare(self, bio_system):
        IFTFlagellaPart.prepare(self, bio_system)
        self.lmb_p = bio_system.getConstantValue("lmb_p")
        self.s_p = bio_system.getConstantValue("s_p")
        self.phi = bio_system.getConstantValue("_phi")

    def process(self, y, t):
        print_log("IFTMovePlus")
        yy = list(y)
        self.ift_wrapper.set_y(yy)
        self.ift_wrapper["X"] += self.s_p
        return yy

    def get_rate(self, y, t):
        self.ift_wrapper.set_y(y)
        flagella_wrapper = self.get_ift_flagella(y)
        if(self.ift_wrapper["X"] <= self.phi * flagella_wrapper["U"] - self.s_p and self.ift_wrapper["D"] == +1):
            return float(self.lmb_p)
        return 0

class IFTMoveMinus(IFTFlagellaPart):
    def __init__(self, part_name, ift_wrapper):
        IFTFlagellaPart.__init__(self, part_name, ift_wrapper)

    def prepare(self, bio_system):
        IFTFlagellaPart.prepare(self, bio_system)
        self.lmb_m = bio_system.getConstantValue("lmb_m")
        self.s_m = bio_system.getConstantValue("s_m")
        self.phi = bio_system.getConstantValue("_phi")

    def process(self, y, t):
        print_log("IFTMoveMinus")
        yy = list(y)
        self.ift_wrapper.set_y(yy)
        self.ift_wrapper["X"] -= self.s_m
        return yy

    def get_rate(self, y, t):
        self.ift_wrapper.set_y(y)
        flagella_wrapper = self.get_ift_flagella(y)
        if(self.ift_wrapper["X"] >= self.s_m and self.ift_wrapper["D"] == -1):
                return float(self.lmb_m)
        return 0


class IFTDetachToZoneB(IFTFlagellaPart):
    def __init__(self, part_name, ift_name):
        IFTFlagellaPart.__init__(self, part_name, ift_name)

    def prepare(self, bio_system):
        IFTFlagellaPart.prepare(self, bio_system)
        self.lmb_2 = bio_system.getConstantValue("lmb_2")
        self.s_p = bio_system.getConstantValue("s_p")
        self.phi = bio_system.getConstantValue("_phi")

    def process(self, y, t):
        print_log("IFTDetachToZoneB")
        yy = list(y)
        self.ift_wrapper.set_y(yy)
        flagella_wrapper = self.get_ift_flagella(yy)
        self.ift_wrapper["X"] = -1
        self.ift_wrapper["D"] = +1
        flagella_wrapper["U"] += self.ift_wrapper["C"]
        self.ift_wrapper["C"] = 0
        return yy

    def get_rate(self, y, t):
        self.ift_wrapper.set_y(y)
        flagella_wrapper = self.get_ift_flagella(y)
        if(self.ift_wrapper["X"] > self.phi * flagella_wrapper["U"] - self.s_p and self.phi * self.ift_wrapper["D"] == +1):
            return float(self.lmb_2)
        return 0


class IFTZoneBAttachPart(IFTFlagellaPart):
    def __init__(self, part_name, ift_wrapper):
        IFTFlagellaPart.__init__(self, part_name, ift_wrapper)

    def prepare(self, bio_system):
        IFTFlagellaPart.prepare(self, bio_system)
        self.lmb_3 = bio_system.getConstantValue("lmb_3")
        self.c_max = bio_system.getConstantValue("C_max")
        self.phi = bio_system.getConstantValue("_phi")

    def process(self, y, t):
        print_log("IFTZoneBAttachPart")
        yy = list(y)
        self.ift_wrapper.set_y(yy)
        flagella_wrapper = self.get_ift_flagella(yy)
        self.ift_wrapper["X"] = self.phi * flagella_wrapper["U"]
        self.ift_wrapper["D"] = -1
        self.ift_wrapper["C"] = 0
        return yy

    def get_rate(self, y, t):
        self.ift_wrapper.set_y(y)
        if(self.ift_wrapper["X"] == -1 and self.ift_wrapper["D"] == +1):
            return float(self.lmb_3)
        return 0

class IFTDetachToZoneA(IFTFlagellaPart):
    def __init__(self, part_name, ift_wrapper):
        IFTFlagellaPart.__init__(self, part_name, ift_wrapper)

    def prepare(self, bio_system):
        IFTFlagellaPart.prepare(self, bio_system)
        self.lmb_4 = bio_system.getConstantValue("lmb_4")
        self.s_m = bio_system.getConstantValue("s_m")
        self.phi = bio_system.getConstantValue("_phi")
        self.a_zone_index = bio_system.compositorIndex("A")

    def process(self, y, t):
        print_log("IFTDetachToZoneA")
        yy = list(y)
        self.ift_wrapper.set_y(yy)
        self.ift_wrapper["X"] = -1
        self.ift_wrapper["D"] = -1
        return yy

    def get_rate(self, y, t):
        self.ift_wrapper.set_y(y)
        print_log(self.s_m)
        print_log(self.ift_wrapper["X"] <= self.s_m)
        if(self.ift_wrapper["X"] <= self.s_m and self.ift_wrapper["D"] == -1):
            return float(self.lmb_4)
        return 0


class BioSystem_StochasticCellFlagella(pbfs.BioSystem):

    def __init__(self, init_sys):
        # Create a BioSystem to simulate.
        pbfs.BioSystem.__init__(self)
        init_sys(self)
        self.flagella_setup()

    def flagella_setup(self):
        sys = self
            # A zone element count
        c = sys.addCompositor("A", 0)

        # IFT position compositors vector
        self.ift_x = {}
        self.ift_d = {}
        self.ift_c = {}
        self.ift_z = {}
        self.ift_wrappers = []
        for i in range(0, sys.getConstantValue("_N")):
            ift_name = "ift_" + str(i)
            x_name = ift_name + "_X"
            x = sys.addCompositor(x_name, -1) # Initial value for all IFT state
            self.ift_x[x_name] = x

            d_name = ift_name + "_D"
            d = sys.addCompositor(d_name, -1) # Initial value for all IFT state
            self.ift_d[d_name] = d

            c_name = ift_name + "_C"
            c = sys.addCompositor(c_name, 0) # Initial value for all IFT state
            self.ift_c[c_name] = c

            z_name = ift_name + "_Z"
            z = sys.addCompositor(z_name, 0) # Initial value for all IFT state
            self.ift_z[z_name] = z

            self.ift_wrappers.append(SubCompositorWrapper(self, 'ift', i))

        # Flagella initial state
        self.flagella_L = {}
        self.flagella_B = {}
        self.flagella_wrappers = []
        for j in range(0, sys.getConstantValue("_M")):
            flagella_name = "flagella_" + str(j)
            l_name = flagella_name + "_U"
            l = sys.addCompositor(l_name, 0) # Initial value for all Flagella state
            self.flagella_L[l_name] = l
            self.flagella_wrappers.append(SubCompositorWrapper(self, 'flagella', j))


        IFT_move_parts = {}
        for i in range(0, sys.getConstantValue("_N")):
            ift_compositor_name = "ift_" + str(i)
            ift_compositor_wrapper = SubCompositorWrapper(self, 'ift', i)
            # IFT attach to flagella
            for j in range(0, sys.getConstantValue("_M")):
                flagella_compositor_name = "flagella_" + str(j)
                part_name = ift_compositor_name + "_attach_" + flagella_compositor_name
                attach_part = IFTZoneAAttachPart(part_name, ift_compositor_wrapper, SubCompositorWrapper(self, 'flagella', j))
                sys.addPart(attach_part)

            # IFT move plus
            part_name = ift_compositor_name + "_move_plus" 
            move_plus_part = IFTMovePlus(part_name, ift_compositor_wrapper)
            sys.addPart(move_plus_part)

            #IFT detach to zone B
            part_name = ift_compositor_name + "_detach_to_b" 
            detach_to_b_part = IFTDetachToZoneB(part_name, ift_compositor_wrapper)
            sys.addPart(detach_to_b_part)

            #IFT attach to flagella at zone B
            part_name = ift_compositor_name + "_attach_b" 
            attach_zone_b_part = IFTZoneBAttachPart(part_name, ift_compositor_wrapper)
            sys.addPart(attach_zone_b_part)

            # IFT move plus
            part_name = ift_compositor_name + "_move_minus" 
            move_minus_part = IFTMoveMinus(part_name, ift_compositor_wrapper)
            sys.addPart(move_minus_part)

            #IFT detach to zone A
            part_name = ift_compositor_name + "_detach_to_a" 
            detach_to_a_part = IFTDetachToZoneA(part_name, ift_compositor_wrapper)
            sys.addPart(detach_to_a_part)

        for j in range(0, sys.getConstantValue("_M")):
            #IFT detach to zone A
            flagella_compositor_name = "flagella_" + str(j)
            part_name = flagella_compositor_name + "_degradation" 
            flagella_degradation_part = DegradationPart(part_name, SubCompositorWrapper(self, 'flagella', j))
            sys.addPart(flagella_degradation_part)

        return sys