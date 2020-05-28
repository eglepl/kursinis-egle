import sys
sys.path.append('..')

import numpy as np
from os import *
from os.path import isfile, isdir, join
import matplotlib.pyplot as plt
from utils import *

file1 = "../output/pandas-example_cell_flagellas_M_N_ode_5.py(3,0)_20200528T121205_35332"
file2 = "../output/pandas-example_cell_flagellas_M_N_stochastic_5.py(4,40)_20200528T020150_38499"

pairs = [
    ("L0", "flagella_0_U"),
    ("L1", "flagella_1_U"),
    ("L2", "flagella_2_U"),
    ("A", "A"),
]

(f1const, f1data) = read_csv(file1)
(f2const, f2data) = read_csv(file2)

# print(f2data)

T = np.linspace(0, 2000, 2000)

correlations = []
for p in pairs:
    (f1key, f2key) = p
    f1_vals = np.interp(T, f1data["time"].tolist(), f1data[f1key].tolist())
    f2_vals = np.interp(T, f2data["time"].tolist(), f2data[f2key].tolist())

    correlations.append((f1_vals, f2_vals))

plt.figure()
# Create a 5% (0.05) and 10% (0.1) padding in the
# x and y directions respectively.
plt.margins(0.05, 0.1)
plt.grid(True)

#plt.plot(t, t_meandata[:, 1], label="A zone mean")
#plt.plot(t, t_stddata[:, 1], label="A zone std")

for (x, y) in correlations:
    plt.scatter(x, y, label="", alpha=0.3)

np.polyfit(correlations[0][0], correlations[0][1], 1)

fig, axs = plt.subplots(1, 4, figsize=(8, 2))
axs[0].scatter(correlations[0][0], correlations[0][1], label='Flagella #1 length')
axs[1].scatter(correlations[1][0], correlations[1][1], label='Flagella #2 length')
axs[2].scatter(correlations[2][0], correlations[2][1], label='Flagella #3 length')
axs[3].scatter(correlations[3][0], correlations[3][1], label='Zone A elements')
fig.suptitle('ODE and Stochastics scatter plot, when zone A is empty')

axs[0].legend()
axs[0].set_xlabel('ODE, units of length')
axs[0].set_ylabel('Stochastics, units of length')

axs[1].legend()
axs[1].set_xlabel('ODE, units of length')
axs[1].set_ylabel('Stochastics, units of length')

axs[2].legend()
axs[2].set_xlabel('ODE, units of length')
axs[2].set_ylabel('Stochastics, units of length')

axs[3].legend()
axs[3].set_xlabel('ODE, elements')
axs[3].set_ylabel('Stochastics, elements')

plt.show()
