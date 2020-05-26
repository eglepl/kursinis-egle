# -*- coding: utf-8 -*-

from sympy import *
from Part import Part

## Representation of a chemical reaction.
#
#  A Part is a process, changing the values of compositors according
#  to some rate laws.
#
#  @author Eglė Plėštytė
#  @date 2017-05-10

class SimplePart(Part):


    ## The constructor
    #  @param self The object pointer.
    #  @param name Name assigned to a StochasticPart.
    #  @param compositors Compositors involving part.
    #  @param change_vector Stochastic chemical reaction molecule changes
    # vector.
    #  @param rate Part rate.
    def __init__(self, name, compositors, change_by, rate):
        ## Name assigned to a Part.
        self.name = name

        self.compositors = compositors
        self.compositors_sys_indexes = None

        self.change_by = change_by
        self.change_byf = None

        self.rate = rate
        self.ratef = None


    def prepare(self, bio_system):
        const_substitutions = bio_system.constant_substitutions

        rate_exp = sympify(self.rate).subs(const_substitutions)
        self.ratef = lambdify(bio_system.symbols, rate_exp)

        self.change_byf = [lambdify(bio_system.symbols,
                        sympify(change).subs(const_substitutions))
                        for change in self.change_by]

        self.compositors_sys_indexes = []
        for c in self.compositors:
            self.compositors_sys_indexes.append(bio_system.compositors.index(c))


    ## Modify the state when this part stochastic event happens.
    #  @param self The object pointer.
    #  @param y current values of the BioSystem Compositors
    #  @param constants dictionary of constants
    #  @param bio_system the BioSystem that is simulated
    #  @return list of updated values of BioSystem compositors
    def process(self, y, t):
        yy = list(y)
        ty = [t] + yy
        for i in range(0, len(self.compositors)):
            compositor_sys_index = self.compositors_sys_indexes[i]
            compositor_change_byf = self.change_byf[i]

            yy[compositor_sys_index] += compositor_change_byf(*ty)
        return yy


    ## Get the stochastic rate of this part using current state of BioSystem
    #  @param self The object pointer.
    #  @param y current values of the BioSystem Compositors
    #  @param constants dictionary of constants
    #  @param bio_system the BioSystem that is simulated
    #  @return stochastic rate of current part
    def get_rate(self, y, t):
        ty = [t] + list(y)
        rate = self.ratef(*ty)
        return rate
