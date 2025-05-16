import numpy as np
import matplotlib.pyplot as plt
import read

t = read.read_0d('t.dac')
x = read.read_0d('x.dac')
vx = read.read_1d('vx.dac')
rho = read.read_1d('rho.dac')
p = read.read_1d('p.dac')

#print(max(vx[0]))


# for graphing
fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(111)
#ax.set_aspect('equal')
plt.rcParams['font.size']=15
plt.tick_params(labelsize=10)


i = 6


plt.plot(x, vx[i],label='vx')
plt.plot(x, p[i],label='p')
plt.plot(x,rho[i],label='rho')
plt.title(str(t[i])[:4])

plt.legend()
plt.show()
