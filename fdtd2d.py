from ctypes import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

file1 = CDLL('./fdtd2d.dll')

file1.fdtd2d.argtypes = [c_double, c_double,c_double,c_double,c_double,c_double,c_float,c_float,c_double,c_double,c_double,c_double,c_double,c_double,c_double,c_double]
file1.fdtd2d.restype = c_uint

b = 0.00000005
l = b

run = file1.fdtd2d(300e-9,300e-9,2e-9,2e-9,10,1e6,25e-9,l,b,-2,-2 , 1, 1e8, 1, 1,80e-9)

for i in range(1,1000):
    data = pd.read_csv("./data/Edata" + str(i) + ".txt")
    plt.imshow(data, cmap = 'RdBu', interpolation='nearest')
    plt.show()
    plt.pause(0.5)