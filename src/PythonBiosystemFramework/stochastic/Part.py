# -*- coding: utf-8 -*-

## Representation of a stochastic process event.
#
#  A Part is a process, changing the values of compositors according
#  to some process and with rates.
#
#  @author Eglė Plėštytė
#  @date 2017-05-10

class Part:

    ## The constructor
    #  @param self The object pointer.
    #  @param name Name assigned to a StochasticPart.
    def __init__(self, name):
        ## Name assigned to a Part.
        self.name = name


    ## Prepare the part for simulations, method is called before BioSystem
    # simulation loop. Every preparation must be done here
    # when it is known what BioSystem will simulate the part
    #   @param bio_system Stochastic BioSystem that will use process and
    # get_rate methods for some time from now.
    #   @return None
    def prepare(self, bio_system):
        pass


    ## Modify the state when this part stochastic event happens.
    #  @param self The object pointer.
    #  @param y current values of the BioSystem Compositors
    #  @return list of updated values of BioSystem compositors, might be a list
    # of sympy expressions (evaluated with BioSystem Compositor, Constants
    # values results updated value).
    def process(self, y):
        raise NotImplementedError


    ## Get the stochastic rate of this part using current state of BioSystem
    #  @param self The object pointer.
    #  @param y current values of the BioSystem Compositors
    #  @return stochastic rate of current part, might be sympy expression (
    # evaluated with BioSystem values results rate).
    def get_rate(self, y):
        raise NotImplementedError
