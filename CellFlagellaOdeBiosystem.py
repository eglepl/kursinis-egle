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

class BioSystem_OdeCellFlagella:

    def __init__(self, init_sys):
        # Create a BioSystem to simulate.
        self.sys = pbfo.BioSystem()
        init_sys(self.sys)
        self.setup()

    def get_sys(self):
        return self.sys

    def setup(self):
        sys = self.sys

        # A zone element count
        A = sys.addCompositor('A', sys.getConstantValue('A_zone_init'))

        # Flagella initial state
        L = list()
        for j in range(0, sys.getConstantValue('_M')):
            c_name = 'L' + str(j)
            c = sys.addCompositor(c_name, 0) # Initial value for all Flagella state
            L.append(c)

        B = list()
        for j in range(0, sys.getConstantValue('_M')):
            c_name = 'B' + str(j)
            c = sys.addCompositor(c_name, 0) # Initial value for all Flagella state
            B.append(c)

        for j in range(0, sys.getConstantValue('_M')):
            reaction  = pbfo.Part(
            'L increase',
            [L[j], A],
            [pbfo.Rate('_beta * (3 / ( (1 / (_lmb_1 / (1 + _alpha * L'+str(j)+'))) + (1/_lmb_p) + (1/_lmb_2) )) * A'), 
            pbfo.Rate('-_beta * (3 / ( (1 / (_lmb_1 / (1 + _alpha * L'+str(j)+'))) + (1/_lmb_p) + (1/_lmb_2) )) * A')])
            sys.addPart(reaction)

        for j in range(0, sys.getConstantValue('_M')):
            reaction  = pbfo.Part(
            'Tip decompose',
            [L[j], B[j]],
            [pbfo.Rate('-_mu'), pbfo.Rate('_mu')])
            sys.addPart(reaction)

        for j in range(0, sys.getConstantValue('_M')):
            reaction  = pbfo.Part(
            'B decrease, A increase',
            [B[j], A],
            [pbfo.Rate('-1 * _beta * (3/ ( (1/_lmb_3) + (1/_lmb_m) + (1/_lmb_4) )) * B' + str(j)), 
            pbfo.Rate('_beta * (3/ ( (1/_lmb_3) + (1/_lmb_m) + (1/_lmb_4) )) * B' + str(j))])
            sys.addPart(reaction)

        return sys