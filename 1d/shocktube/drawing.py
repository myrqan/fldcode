import numpy as np
import matplotlib.pyplot as plt
import read

t = read.read_0d('t.dac')
x = read.read_0d('x.dac')
vx = read.read_1d('vx.dac')
rho = read.read_1d('rho.dac')
p = read.read_1d('p.dac')

# for graphing
fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(111)
ax.set_aspect('equal')
plt.rcParams['font.size']=15
plt.tick_params(labelsize=10)




plt.plot(x, vx[3],label='vx')
plt.plot(x, p[3],label='p')
plt.plot(x,rho[3],label='rho')
plt.title(str(t[3])[:4])

plt.legend()
plt.show()
