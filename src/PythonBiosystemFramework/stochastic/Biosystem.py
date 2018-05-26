# -*- coding: utf-8 -*-

import numpy as np
from sympy import *
from BiosystemBase import *
from Compositor import Compositor
from Const import Const


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

class StochasticBioSystem(BioSystemBase):

    ## The constructor
    #  @param self The object pointer.
    #  @param compositors Compositors involving chemical reactions.
    #  @param constants Constants of a Biosystem.
    #  @param symbols A Symbol with a @p name representing a substance.
    def __init__(self):
        ## A seed for random number generator for reproduction of results
        self.seed = None

        ## A seed for random number generator for reproduction of results
        self.part_rate_functions = None

    ## Simulation of the Biosystem stochasticaly.
    #  @param self The object pointer.
    #  @param tspan Time interval to simulate, for example [t0, t1],
    #   if contains more than 2 elements it will be used as time points
    #   to calculate values.
    #  @param sample_count (optional) number of points to sample solutioin,
    #   for example 100, default value is 100
    #  @return Tuple (T, Y), where T - time point list, Y - matrix consisting
    #  of Compositor values at a time points.
    def _run(self, sample_at):
        y0 = []
        for c in self.compositors:
            y0.append(c.value)

        self.determine_a_mu_functions()

        t = sample_at[0]
        result_y = np.zeros((0, len(self.compositors)))
        sample_i = 0


        while sample_i < len(sample_at):
            if sample_at[sample_i] <= t:
                result_y = np.concatenate([result_y, np.array([y0])])
                sample_i = sample_i + 1

            else:
                gallipie_result = self.gallipie(y0)
                if gallipie_result:
                    (s_dt, part_idx) = gallipie_result
                    curr_part = self.parts[part_idx]
                    y0 = np.add(y0, curr_part.rates)
                    t = t + s_dt
                else:
                    t = sample_at[-1]

        self.result = (sample_at, result_y)
        return self.result

    def gallipie(self, y):
        a_mu = self.determine_a_mu_values(y)
        a = sum(a_mu)
        if not a:
            return None

        s_dt = np.random.exponential(1.0 / a)
        part_idx = np.random.choice(len(a_mu), 1, p=[mu/a for mu in
                                                     a_mu]).item(0)
        return s_dt, part_idx

    def determine_a_mu_functions(self):
        a_mu = []
        for p in self.parts:
            ## All Constant symbols.
            constant_syms = list(map((lambda c: c.sym), self.constants))
            ## All Constant values.
            constant_vals = list(map((lambda c: c.value), self.constants))
            ## Pairs of symbols and its value.
            substitutions = zip(constant_syms, constant_vals)
            prob_coeff_expr = sympify(p.prob_coeff)
            prob_coeff_expr = prob_coeff_expr.subs(substitutions)
            a = lambdify(self.symbols, prob_coeff_expr)
            a_mu.append(a)
        self.a_mu_functions = a_mu

    def determine_a_mu_values(self, y):
        a_mu_values = []
        yy = list(y)
        yy.insert(0, 0)
        for a_mu_fn in self.a_mu_functions:
            val = a_mu_fn(*yy)
            a_mu_values.append(val)
        return a_mu_values