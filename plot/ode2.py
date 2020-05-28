import sys
sys.path.append('..')

import matplotlib.pyplot as plt
from utils import *

(const, data) = read_csv("../output/pandas-example_cell_flagellas_M_N_ode_2.py(2,0)_20200528T115516_96537")


# print(const)
# print(data['time'].tolist())

T = data['time'].tolist()
FU = data['L0'].tolist()
FU2 = data['L1'].tolist()
FL = [u for u in FU ]
FL2 = [u for u in FU2 ]
A = data['A'].tolist()

# print(FL)

# plt.figure(num=__file__)
# Create a 5% (0.05) and 10% (0.1) padding in the
# x and y directions respectively.
# plt.margins(0.05, 0.1)
# plt.grid(True)

# plt.plot(T, A, label="A zone")
# plt.plot(T, FL, label="Flagellar length")


fig, axs = plt.subplots(1, 2, figsize=(10, 5))
axs[0].plot(T, A, label="A zone")
axs[1].plot(T, FL, label="Flagellar #1 length")
axs[1].plot(T, FL2, label="Flagellar #2 length")
fig.suptitle('Flagella length dynamics when is deficit in zone A of flagella proteins')

axs[0].legend()
axs[0].set_xlabel('Time')
axs[0].set_ylabel('Elements')

axs[1].legend()
axs[1].set_xlabel('Time')
axs[1].set_ylabel('Units of length')

plt.show()