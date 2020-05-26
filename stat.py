import numpy as np
from os import *
from os.path import isfile, isdir, join
import matplotlib.pyplot as plt

analyze = [
    {"file": 'example_cell_flagellas_M_N_stochastic_2', "cols": [0, 1, 82, 84]},
    {"file": 'example_cell_flagellas_M_N_stochastic_3', "cols": [0, 1, 82, 84]},
    {"file": 'example_cell_flagellas_M_N_stochastic_4', "cols": [0, 1, 42]},
    {"file": 'example_cell_flagellas_M_N_stochastic_5', "cols": [0, 1, 162, 164, 166]},
]

output_dir = './output'

files = []
for (dirpath, dirnames, filenames) in walk('output', True, None):
    files.extend(filenames)
    break

t = np.linspace(0, 1000, 1000)

for experiment in analyze:
    print("experiment: " + experiment['file'])
    matching = filter(lambda fname: fname.startswith(experiment["file"]), files)
    ex_data = []
    for m in matching:
        data = np.genfromtxt('output/' + m, delimiter=';')
        ex_data.append(data)

    nex_data = []
    for data_table in ex_data:
        table = None
        for col_idx in range(1, data_table.shape[1]):
            coly = np.interp(t, data_table[:,0], data_table[:,col_idx])
            if table is None:
                table = coly
            else:
                table = np.vstack((table, coly))
        normalized_data_table = np.transpose(table)
        nex_data.append(normalized_data_table)

    data3darray = np.dstack(nex_data)

    meandata = np.mean(data3darray, axis=2)
    stddata = np.std(data3darray, axis=2)

    t_meandata = np.concatenate((np.vstack(t), meandata), axis=1)
    t_stddata = np.concatenate((np.vstack(t), stddata), axis=1)


    np.savetxt(join("stat", experiment["file"] + ".mean.stat"), t_meandata[:, experiment["cols"]], delimiter=';', fmt='%10.5f')
    np.savetxt(join("stat", experiment["file"] + ".std.stat"), t_stddata[:, experiment["cols"]], delimiter=';', fmt='%10.5f')


    plt.figure(num=experiment["file"])
    # Create a 5% (0.05) and 10% (0.1) padding in the
    # x and y directions respectively.
    plt.margins(0.05, 0.1)
    plt.grid(True)

    #plt.plot(t, t_meandata[:, 1], label="A zone mean")
    #plt.plot(t, t_stddata[:, 1], label="A zone std")

    plt.plot(t, t_meandata[:, 2], label="Flagella #1 length mean")
    plt.plot(t, t_stddata[:, 2], label="Flagella #1 length std")

    plt.legend()
    plt.xlabel('Time')
    plt.ylabel('Elements')
    plt.show()