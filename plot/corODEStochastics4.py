import sys
sys.path.append('..')

import numpy as np
from os import *
from os.path import isfile, isdir, join
import matplotlib.pyplot as plt
from utils import *

file1 = "../output/pandas-example_cell_flagellas_M_N_ode_4.py(1,0)_20200528T120716_1761"
file2 = "../output/pandas-example_cell_flagellas_M_N_stochastic_4.py(1,10)_20200528T012436_15310"

pairs = [
    ("L0", "flagella_0_U"),
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

fig, axs = plt.subplots(1, 2, figsize=(10, 5))
axs[0].scatter(correlations[0][0], correlations[0][1], label='Flagella #1 length')
axs[1].scatter(correlations[1][0], correlations[1][1], label='Zone A elements')
fig.suptitle('ODE and Stochastics scatter plot. Initial flagella growth dynamics.')

axs[0].legend()
axs[0].set_xlabel('ODE, units of length')
axs[0].set_ylabel('Stochastics, units of length')

axs[1].legend()
axs[1].set_xlabel('ODE, elements')
axs[1].set_ylabel('Stochastics, elements')

plt.show()
