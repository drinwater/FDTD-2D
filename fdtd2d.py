from ctypes import *
from ssl import ALERT_DESCRIPTION_CLOSE_NOTIFY
from timeit import repeat
from tkinter import Canvas
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

file1 = CDLL('./fdtd2d.dll')

file1.fdtd2d.argtypes = [c_double, c_double,c_double,c_double,c_double,c_double,c_float,c_float,c_double,c_double,c_double,c_double,c_double,c_double,c_double,c_double]
file1.fdtd2d.restype = c_uint

b = 0.00000005
l = b

fig, ax = plt.subplots()
i=1
print(str(i)+'\n')
#ax.imshow(data, cmap = 'RdBu', interpolation='nearest', vmin=-0.3, vmax=0.3)
run = file1.fdtd2d(300e-9,300e-9,2e-9,2e-9,10,1e6,25e-9,50e-9,50e-9,-2,-2, 1, 1e8, 1, 1,80e-9)

def plotfunc(i):
    for i in range(1,1000):
        print(str(i)+'\n')
        data = np.loadtxt("./data/Edata" + str(i) + ".txt")
        ax.imshow(data, cmap = 'RdBu', interpolation='nearest', vmin=-0.3, vmax=0.3)


anim = FuncAnimation(fig, plotfunc, interval=500, frames=10000)
 
plt.draw()
plt.show()