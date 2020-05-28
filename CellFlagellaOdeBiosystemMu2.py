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
        self.Bk = []
        self.Ck = []
        init_sys(self)
        # A zone element count
        self.Ae = self.addCompositor('Ae', 0)
        # A zone kinesin concentration
        self.Ak = self.addCompositor('Ak', 0)

    def add_flagella(self, idx, l, bk):
        # flagella length
        c_name = 'L' + str(idx)
        c = self.addCompositor(c_name, l) # Initial value for all Flagella state
        self.L.append(c)

        # kinesin count in B zone
        c_name = 'Bk' + str(idx)
        c = self.addCompositor(c_name, bk) # Initial value for all Flagella state
        self.Bk.append(c)


    def setup(self):
        for j in range(0, len(self.L)):
            tip_construct  = pbfo.Part(
                'tip_construct',
                [self.L[j], self.Bk[j], self.Ak],
                [pbfo.Rate('(_beta * Ae) * _N * _omega  / (1 + 2 * L' + str(j) +')'),
                pbfo.Rate('_gamma * _N * _omega * Ak  / (1 + 2 * L' + str(j) +')'),
                pbfo.Rate('-_gamma * _N * _omega * Ak  / (1 + 2 * L' + str(j) +')')
                ])
            self.addPart(tip_construct)

            # delta = length of tip, _eta - area of the tip,  _delta * _eta - volume of the tip
            tip_deconstruct  = pbfo.Part(
                'tip_deconstruct',
                [self.L[j]],
                [pbfo.Rate('-1 * _mu * Bk' + str(j))])
            self.addPart(tip_deconstruct)

            kinesin_diffuse_BA  = pbfo.Part(
                'kinesin_diffuse_BA',
                [self.Bk[j], self.Ak],
                [pbfo.Rate('-_kba * (Bk' + str(j) + ' - Ak) / (1 + L' + str(j) + ")"),
                pbfo.Rate('_kba * (Bk' + str(j) + ' - Ak) / (1 + L' + str(j) + ")")])
            self.addPart(kinesin_diffuse_BA)

        return sys