import numpy as np
import matplotlib.pyplot as plt
import os


xx = np.fromfile('dat/x.dat',dtype=np.float64)
ro = np.fromfile('dat/ro.dat',dtype=np.float64)

## for draw graph
fig = plt.figure(figsize=(7, 5))
plt.rcParams['font.size']=13
plt.rcParams['font.family']='STIXGeneral'
plt.rcParams['mathtext.fontset']='stix'
ax = fig.add_subplot(111)


ax.plot(xx,ro)

plt.show()

