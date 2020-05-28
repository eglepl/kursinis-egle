import sys
import numpy as np
import pandas as pd
from datetime import datetime
from random import *

def print_log(s):
    #print(s)
    None

def save_csv(filename, sys, TY):

    filename = "output/pandas-" + filename + "(" + str(sys.getConstantValue('_M')) + "," + str(sys.getConstantValue('_N')) + ")"
    filename += "_" + datetime.now().strftime("%Y%m%dT%H%M%S")
    filename += "_" + str(randrange(100000))

    data_frame_const = pd.DataFrame(data=[[c.value for c in sys.constants]],    # values
              columns=[c.name for c in sys.constants])  # 1st row as the column names

    data_frame = pd.DataFrame(data=TY,    # values
              columns=['time'] + [c.name for c in sys.compositors])  # 1st row as the column names

    data_frame.to_csv(filename + ".data.csv")
    data_frame_const.to_csv(filename + ".const.csv")

def read_csv(filename):
    data_frame = pd.read_csv(filename + ".data.csv")
    data_frame_const = pd.read_csv(filename + ".const.csv")

    constants = {}
    for i in data_frame_const.columns:
        constants[i] = data_frame_const.iloc[0][i]

    return (constants, data_frame)

