# -*- coding: utf-8 -*-

import numpy as np
import PythonBiosystemFramework.BiosystemBase as b
from PythonBiosystemFramework.stochastic.Compositor import Compositor
import random

## Biological system to simulate stochastic
#
#  In order to analyze a biological system you create a BioSystem
#  object.
#
#  Biological system might need some constants (Const) and compositors
#  (Compositor). Compositors are the total rates of change of substance state
#  variables.
#
#  Then parts (Part) are declared and added to a system. Part is a process in
#  a system (reaction) that affect some state variables by changing their
#  value  according to a rate law (Rate).
#
#  Then a simulation can be started.
#
#  Example:
#
#  @code
#
# from Biosystem import *
# from Part import *
# from Rate import *
# from Pulse import *
# import matplotlib.pyplot as plt
#
# sys = BioSystem()
# sys.addConstant('k', 0.05)
# dAdt = sys.addCompositor('A', 10)
# dBdt = sys.addCompositor('B', 0)
# dEdt = sys.addCompositor('E', 1)
# reaction  = Part(
# 'A + E -k> B + E',
# [dAdt, dBdt, dEdt],
# [Rate('-k * A * E'), Rate('k * A * E'), Rate('0')])
# sys.addPart(reaction)
# T = None
# Y = None
# (T, Y) = sys.run([0, 25])
#
# # Plot the simulation data
# plt.figure()
# plt.plot(T, Y[:, sys.compositorIndex('A')], label="A")
# plt.plot(T, Y[:, sys.compositorIndex('B')], label="B")
# plt.plot(T, Y[:, sys.compositorIndex('E')], label="E")
# plt.legend()
# plt.xlabel('Time')
# plt.ylabel('Concentration')
# plt.show()
#
#  @endcode
#
#  @author Eglė Plėštytė
#  @date 2017-06-15


class BioSystem(b.BioSystemBase):

    ## The constructor
    #  @param self The object pointer.
    #  @param compositors Compositors involving chemical reactions.
    #  @param constants Constants of a Biosystem.
    #  @param symbols A Symbol with a @p name representing a substance.
    def __init__(self):
        b.BioSystemBase.__init__(self)
        ## A seed for random number generator for reproduction of results
        self.seed = None
        self.max_accuracy = False
        self.constant_substitutions = None

    ## Create or add a compositor to the system
    #  @param self The object pointer.
    #  @param compositor_or_name
    #     Compositor object to add if @p init_value == None,
    #     else Compositor name to add with an @p init_value.
    #
    #  @param init_value Initial value of a Compositor or None.
    #  @return The Compositor object that was added.
    def addCompositor(self, compositor_or_name, init_value=None):
        new_compositor = None
        if init_value == None:
            new_compositor = compositor_or_name
        else:
            name = compositor_or_name
            new_compositor = Compositor(name, init_value)

        b.BioSystemBase.addCompositor(self, new_compositor)

        return new_compositor

    ## Simulation of the Biosystem stochasticaly.
    #  @param self The object pointer.
    #  @param sample_at Sample the simulation at provided time points.
    #  @return Tuple (T, Y), where T - time point list, Y - matrix consisting
    #  of Compositor values at a time points.
    def _run(self, sample_at):

        self.setup_constant_substitutions()
        for p in self.parts:
            p.prepare(self)

        y0 = []
        y0_before = None
        for c in self.compositors:
            y0.append(c.value)

        sample_i = 0
        t = sample_at[0]
        t_end = sample_at[-1]

        # result_t = np.empty((0, 1))
        # result_y = np.empty((0, len(y0)))
        result_t = []
        result_y = []

        while True:
            if self.max_accuracy:
                if y0_before:
                    # result_t = np.vstack([result_t, [t]])
                    # result_y = np.vstack((result_y, y0_before))
                    result_t.append(t)
                    result_y.append(y0_before)
                # result_t = np.vstack([result_t, [t]])
                # result_y = np.vstack((result_y, y0))
                result_t.append(t)
                result_y.append(y0)
            else:
                while (sample_i < len(sample_at)) and (t >= sample_at[sample_i]):
                    #result_y = np.vstack((result_y, y0))
                    result_y.append(y0)
                    sample_i += 1

            if t >= t_end:
                break

            gillespie_result = self.gillespie(y0, t)
            if gillespie_result:
                (s_dt, part_idx) = gillespie_result
                part = self.parts[part_idx]
                y0_before = list(y0)
                y0 = part.process(y0, t)
                t = t + s_dt
            else:
                y0_before = list(y0)
                t = sample_at[-1]

        result_y_array = np.array(result_y)

        if self.max_accuracy:
            self.result = (result_t, result_y_array)
        else:
            self.result = (sample_at, result_y_array)

        return self.result

    def gillespie(self, y, t):
        parts_a = self.determine_parts_rates(y, t)
        a = float(sum(parts_a))
        if not a:
            return None

        s_dt = np.random.exponential(1.0 / a)
        #part_idx = np.random.choice(len(parts_a), 1, p=[pa / a for pa in
        #                                                parts_a]).item(0)
        part_idx = self.choice([pa / a for pa in parts_a])

        return s_dt, part_idx

    def determine_parts_rates(self, y, t):
        rates = []
        for p in self.parts:
            part_rate = p.get_rate(y, t)
            rates.append(part_rate)
        return rates

    def setup_constant_substitutions(self):
        ## All Constant symbols.
        constant_syms = list(map((lambda c: c.sym), self.constants))
        ## All Constant values.
        constant_vals = list(map((lambda c: c.value), self.constants))
        ## Pairs of symbols and its value.
        self.constant_substitutions = zip(constant_syms, constant_vals)

    def choice(self, p):
        r = random.random()
        psum = 0
        for i in range(0, len(p)):
            psum += p[i]
            if psum >= r:
                return i
        return len(p) - 1