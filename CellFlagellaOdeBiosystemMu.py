import sys
sys.path.append('./src')

import PythonBiosystemFramework as pbf
import PythonBiosystemFramework.ode as pbfo
from scipy.stats import *
import numpy as np
from datetime import datetime
from random import *
from copy import deepcopy

def print_log(s):
    #print(s)
    None

class BioSystem_OdeCellFlagella(pbfo.BioSystem):

    def __init__(self, init_sys):
        # Create a BioSystem to simulate.
        pbfo.BioSystem.__init__(self)
        self.L = []
        self.B = []
        init_sys(self)

    def add_flagella(self, idx, size, b):
        c_name = 'L' + str(idx)
        c = self.addCompositor(c_name, size) # Initial value for all Flagella state
        self.L.append(c)
        c_name = 'B' + str(idx)
        c = self.addCompositor(c_name, b) # Initial value for all Flagella state
        self.B.append(c)

    def setup(self):
        sys = self

        # A zone element count
        self.A = self.addCompositor('A', sys.getConstantValue('A_zone_init'))

        for j in range(0, len(self.L)):
            reaction  = pbfo.Part(
            'L increase',
            [self.L[j], self.A],
            [pbfo.Rate('_beta * (3 / ( (1 /_lmb_1 ) + (1/_lmb_p) + (1/_lmb_2) + (1 /_lmb_3 ) + (1/_lmb_m) + (1/_lmb_4) )) * A / (1 + _k * 2 * L' + str(j) + ')'), 
            pbfo.Rate('-_beta * (3 / ( (1 /_lmb_1 ) + (1/_lmb_p) + (1/_lmb_2) + (1 /_lmb_3 ) + (1/_lmb_m) + (1/_lmb_4) )) * A / (1 + _k * 2 * L' + str(j) + ')')])
            sys.addPart(reaction)

        for j in range(0, len(self.L)):
            reaction  = pbfo.Part(
            'Tip decompose',
            [self.L[j], self.A],
            [pbfo.Rate('-_mu * L' + str(j)), pbfo.Rate('_mu * L' + str(j))])
            sys.addPart(reaction)

        return sys