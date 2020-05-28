import numpy as np
from os import *
from os.path import isfile, isdir, join
import matplotlib.pyplot as plt

def scatter(compare, T):
    pairs = []
    for experiment in compare:
        ode = experiment[0]
        stoch = experiment[1]
        odeTable = np.genfromtxt(ode["file"], delimiter=';')
        stochTable = np.genfromtxt(stoch["file"], delimiter=';')

        odeValues = np.interp(T, odeTable[:,ode["cols"][0]], odeTable[:,ode["cols"][1]])
        stochValues = np.interp(T, stochTable[:,stoch["cols"][0]], stochTable[:,stoch["cols"][1]])

        pairs.append((experiment[2], odeValues, stochValues))


    plt.figure()
    # Create a 5% (0.05) and 10% (0.1) padding in the
    # x and y directions respectively.
    plt.margins(0.05, 0.1)
    plt.grid(True)

    #plt.plot(t, t_meandata[:, 1], label="A zone mean")
    #plt.plot(t, t_stddata[:, 1], label="A zone std")

    for (lb, o, s) in pairs:
        plt.scatter(o, s, label=lb, alpha=0.3)

    plt.legend()
    plt.xlabel('Ode model')
    plt.ylabel('Stochatic model')
    plt.show()


analyze = [
    [{"file": join('output', 'example_cell_flagellas_M_N_ode_2.py_(2)_ode_0-1500_20200525T124438_32789.csv'), "cols": [0, 5]},
    {"file": join('stat','example_cell_flagellas_M_N_stochastic_2.mean.stat'), "cols": [0, 1]},
    "A Zone"
    ], # time - A zone
  
    [{"file": join('output', 'example_cell_flagellas_M_N_ode_2.py_(2)_ode_0-1500_20200525T124438_32789.csv'), "cols": [0, 1]},
    {"file": join('stat','example_cell_flagellas_M_N_stochastic_2.mean.stat'), "cols": [0, 2]},
    "Flagella #1 length"
    ], # time  - F0
   
    [{"file": join('output', 'example_cell_flagellas_M_N_ode_2.py_(2)_ode_0-1500_20200525T124438_32789.csv'), "cols": [0, 3]},
    {"file": join('stat','example_cell_flagellas_M_N_stochastic_2.mean.stat'), "cols": [0, 3]},
    "Flagella #2 length"
    ], #time - F2
]

analyze2 = [
    [{"file": join('output', 'example_cell_flagellas_M_N_ode_4.py(1)_ode_0-1000_20200525T122523_36661.csv'), "cols": [0, 2]},
    {"file": join('stat','example_cell_flagellas_M_N_stochastic_4.mean.stat'), "cols": [0, 1]},
    "A Zone"
    ], # time - A zone
    
    [{"file": join('output', 'example_cell_flagellas_M_N_ode_4.py(1)_ode_0-1000_20200525T122523_36661.csv'), "cols": [0, 1]},
    {"file": join('stat','example_cell_flagellas_M_N_stochastic_4.mean.stat'), "cols": [0, 2]},
    "Flagella length"
    ], # time  - F0
]


T = np.linspace(0, 500, 500)
scatter(analyze, T)
scatter(analyze2, T)

