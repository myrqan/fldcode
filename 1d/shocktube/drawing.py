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
fig = plt.figure(figsize=(7, 5))
plt.rcParams['font.size']=15
plt.rcParams['font.family']='STIXGeneral'
plt.rcParams['mathtext.fontset']='stix'


i = 6


plt.plot(x, vx[i],label=r'$v_x$')
plt.plot(x, p[i],label=r'$p$')
plt.plot(x,rho[i],label=r'$\rho$')
plt.title(r'$t=$' + (str(t[i])[:4]))

plt.legend(loc='upper left', bbox_to_anchor=(1,1))
plt.tight_layout()
plt.show()
