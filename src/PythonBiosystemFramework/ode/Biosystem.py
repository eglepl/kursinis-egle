# -*- coding: utf-8 -*-

from sympy import *
import numpy as np
from scipy.integrate import *
import PythonBiosystemFramework.BiosystemBase as b
from PythonBiosystemFramework.ode.Compositor import Compositor

## Biological system to simulate
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
        # call parent class constructor
        b.BioSystemBase.__init__(self)
        ## Flag if function determine_rates was called.
        #  If False we need to call determine_rates again.
        self.rates_determined = False

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

    ## Set Constant value by Constant name.
    #  @param self The object pointer.
    #  @param name The name of existing Constant.
    #  @param value New value of the Constant.
    #  @return None.
    def changeConstantValue(self, name, value):
        b.BioSystemBase.changeConstantValue(name, value)
        self.reset_rates()
        return None

    ## Determine rates of all compositors unless already determined.
    #  @param self The object pointer.
    #  @return None.
    def determine_rates(self):
        if not self.rates_determined:
            for i in self.parts:
                p = i
                for k in range(0, len(p.compositors)):
                    p.compositors[k].addRate(p.rates[k])
            for k in self.compositors:
                ## All Constant symbols.
                constant_syms = list( map( (lambda c: c.sym), self.constants))
                ## All Constant values.
                constant_vals = list( map( (lambda c: c.value), self.constants))
                ## Pairs of symbols and its value.
                substitutions = zip(constant_syms, constant_vals)
                ## Convert Compositor rate string/formula to sympy expression
                #  and substitute all constants with its values.
                #print k.rate
                #print substitutions
                symexpr = sympify(k.rate).subs(substitutions)
                #print(symexpr)
                ## Convert sympy expression to python function with arguments:
                #  time, compositor1_name, compositor2_name, ...
                k.ratef = lambdify(self.symbols, symexpr)
            self.rates_determined = True
        return None

    ## Reset all Compositor rates to '0'.
    #  @param self The object pointer.
    #  @return None.
    def reset_rates(self):
        self.rates_determined = False
        for i in self.compositors:
            i.rate = '0'
            i.ratef = None
        return None

    ## Simulate the ODE Biosystem.
    #  @param self The object pointer.
    #  @param tspan Time interval to simulate, for example [t0, t1],
    #   if contains more than 2 elements it will be used as time points
    #   to calculate values.
    #  @param sample_at sample the result at time points
    #  @return Tuple (T, Y), where T - time point list, Y - matrix consisting
    #  of Compositor values at a time points.
    def _run(self, sample_at):

        self.determine_rates()
        y0 = []
        for c in self.compositors:
            y0.append(c.value)

        y = odeint(self.sys_ode, y0, sample_at)

        self.result = (sample_at, y)

        return self.result

    ## Ordinary diferential equatation of the system.
    #  @param self The object pointer.
    #  @param y System compositors values.
    #  @param t Time point.
    #  @return change of Compositor values.
    def sys_ode(self, y, t):
        dy = np.zeros(len(y))

        for i in range(0, len(self.compositors)):
            k = self.compositors[i]
            args = [t] + list(y)  # Add time component to y
            dy[i] = k.ratef(*args)  # list as arguments for function
        return dy

    ## Run simulation given pulse list. Last pulse is not simulated.
    #  Each pulse defines time, Compositor, Compositor value to set at
    #  provided time.
    #
    #  When the Pulse time comes Pulse defined Compositor is set to value
    #  provided. Normal simulation is carried on till next Pulse. While
    #  there is next Pulse - action repeats.
    #
    #  First Pulse should start at t = 0 time.
    #  Each next Pulse time should be greater than previous.
    #  Last Pulse is not simulated - it is a stop time.
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
    #
    # pulses = []
    #
    # # initial condition
    # pulses.append(Pulse(0, 'A', 10))
    #
    # # spike in some A
    # pulses.append(Pulse(100, 'A', 20))
    #
    # # spike in a bit less A
    # pulses.append(Pulse(150, 'A', 5))
    #
    # # spike in more A again
    # pulses.append(Pulse(250, 'A', 10))
    #
    # # stop the simulation at time 500 with this empty string as the state
    # # variable parameter
    # pulses.append(Pulse(500, '', 0))
    #
    # # Run pulsed simulation
    # (T, Y) = sys.run_pulses(pulses)
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
    #  @param self The object pointer.
    #  @param pulse_series List of pulse objects.
    #  @return Tuple (T, Y), where T - time point list, Y - matrix consisting
    #  of Compositor values at a time points.

    # def run_pulses(self, pulse_series, sample_count=1000):
    #     num_pulses = len(pulse_series)
    #
    #     T = []
    #     Y = []
    #
    #     prev_start = 0
    #     prev_end = 0
    #
    #     for i in range(0, num_pulses - 1):
    #         pulse = pulse_series[i]
    #         if pulse.compositor_name:
    #             self.compositors[self.map_compositors[pulse.compositor_name]].value = pulse.value
    #
    #         sim_length = pulse_series[i+1].time - pulse.time
    #         tspan = [prev_end, prev_end + sim_length]
    #         (T_sim, Y_sim) = self.run(tspan, sample_count)
    #
    #
    #         T = T + list(T_sim[2:])
    #         if(len(Y) > 0):
    #             Y = np.concatenate((Y, Y_sim[2:]))
    #         else:
    #             Y = Y_sim[2:]
    #
    #         prev_start = pulse.time
    #         prev_end = prev_start + sim_length
    #
    #         for j in range(0, len(self.compositors)):
    #             self.compositors[j].value = Y_sim[-1][j]
    #
    #     self.reset_state_variables()
    #     return (T, Y)
