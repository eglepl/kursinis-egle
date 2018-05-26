# -*- coding: utf-8 -*-

## Representation of a chemical reaction.
#
#  A Part is a process, changing the values of compositors according
#  to some rate laws.
#
#  @author Eglė Plėštytė
#  @date 2017-05-10

class Part:

    ## The constructor
    #  @param self The object pointer.
    #  @param name Name assigned to a StochasticPart.
    #  @param compositors involving chemical reactions.
    #  @param stochastic chemical reaction molecule changes vector (rates).
    #  @param stochastic chemical reaction rate of the reaction (c_mu).
    #  @param number of distinct molecular reactant combinations of the
    # reaction (h_mu)
    def __init__(self, name, compositors):
        ## Name assigned to a Part.
        self.name = name

    ## Modify the state when this part stochastic event happens.
    #  @param self The object pointer.
    #  @param y current values of the BioSystem Compositors
    #  @param constants dictionary of constants
    #  @param bio_system the BioSystem that is simulated
    #  @return list of updated values of BioSystem compositors, might be a list
    # of sympy expressions (evaluated with BioSystem Compositor, Constants
    # values results updated value).
    def process(self, y, constants, bio_system):
        raise NotImplementedError

    ## Get the stochastic rate of this part using current state of BioSystem
    #  @param self The object pointer.
    #  @param y current values of the BioSystem Compositors
    #  @param constants dictionary of constants
    #  @param bio_system the BioSystem that is simulated
    #  @return stochastic rate of current part, might be sympy expression (
    # evaluated with BioSystem values results rate).
    def get_rate(self, y, constants, bio_system):
        raise NotImplementedError
