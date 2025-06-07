import numpy as np

import matplotlib.pyplot as plt
import read


ixjx = read.read_grid2d('ixjx.dac')
ix = ixjx[0]
jx = ixjx[1]
t = read.read_t('t.dac')
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

