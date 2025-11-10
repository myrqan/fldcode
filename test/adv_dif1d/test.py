import numpy as np
import matplotlib.pyplot as plt

## for draw graph
fig = plt.figure(figsize=(7, 5))
plt.rcParams['font.size']=13
plt.rcParams['font.family']='STIXGeneral'
plt.rcParams['mathtext.fontset']='stix'
ax = fig.add_subplot(111)

ax.set_xlim(0,1)



x = np.linspace(0,1,1001)
y = np.exp(-((x-0.5)/0.05)**2)

plt.plot(x,y)
plt.show()
