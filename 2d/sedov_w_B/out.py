import numpy as np

import matplotlib.pyplot as plt
import read



ixjx = read.read_grid2d('ixjx.dat')
ix = ixjx[0]
jx = ixjx[1]

x = read.read_x2d('x.dat',ix,jx)
z = read.read_x2d('z.dat',ix,jx)


print(x)
print(z)

pr= read.read_phys2d('p.dat',ix,jx)
bz = read.read_phys2d('bz.dat',ix,jx)
etot = read.read_phys2d('etot.dat',ix,jx)

## for draw graph
fig = plt.figure(figsize=(7, 5))
plt.rcParams['font.size']=13
plt.rcParams['font.family']='STIXGeneral'
plt.rcParams['mathtext.fontset']='stix'
ax = fig.add_subplot(111)

#ax.set_xlim(0,0.05)

plt.contourf(x,z,etot[0])
plt.show()

exit()
print(x)


exit()

ixjx = read.read_grid2d('ixjx.dac')
ix = ixjx[0]
jx = ixjx[1]
T = read.read_t('t.dac')
x = read.read_x2d('x.dac',ix,jx)
z = read.read_x2d('z.dac',ix,jx)
ro= read.read_phys2d('rho.dac',ix,jx)
vx= read.read_phys2d('vx.dac',ix,jx)
vz= read.read_phys2d('vz.dac',ix,jx)
pr= read.read_phys2d('p.dac',ix,jx)

nd = np.size(t)

nd = 169 - 1
print("nd=", nd + 1)
print(t[nd])

print(x[25:35,25:35])
print(z[25:35,25:35])
print(vx[nd,25:35,25:35])

