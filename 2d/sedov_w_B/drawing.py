import numpy as np
import matplotlib.pyplot as plt
import read
from scipy.io import FortranFile

ixjx = read.read_grid2d('ixjx.dat')
ix = ixjx[0]
jx = ixjx[1]
t = read.read_t('t.dat')
x = read.read_x2d('x.dat',ix,jx)
z = read.read_x2d('z.dat',ix,jx)
ro= read.read_phys2d('ro.dat',ix,jx)
vx= read.read_phys2d('vx.dat',ix,jx)
vy= read.read_phys2d('vy.dat',ix,jx)
vz= read.read_phys2d('vz.dat',ix,jx)
bx= read.read_phys2d('bx.dat',ix,jx)
by= read.read_phys2d('by.dat',ix,jx)
bz= read.read_phys2d('bz.dat',ix,jx)
pr= read.read_phys2d('p.dat',ix,jx)



X, Z = np.meshgrid(x[:,0],z[0,:])


print(t)
nd = np.size(t)
#print(nd)


## for draw graph
fig = plt.figure(figsize=(10, 8))
plt.rcParams['font.size']=10
plt.rcParams['font.family']='STIXGeneral'
plt.rcParams['mathtext.fontset']='stix'

ax = fig.add_subplot(111)
ax.set_aspect('equal')



n = 0
for n in range(nd):
    ax.set_xlim(0,1)
    ax.set_ylim(-0.1,1)
    im1 = ax.pcolormesh(x,z,ro[n],cmap='plasma',vmin=0,vmax=3)
    im2 = ax.streamplot(X,Z,bx[n],bz[n],density=0.75,color='w')
    cbar1 = fig.colorbar(im1,ax=ax)

    savename = 'fig/'+str(n).zfill(3)+'.png'
    plt.savefig(savename,dpi=300)
    cbar1.remove()
    ax.clear()


exit()
