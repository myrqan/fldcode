import numpy as np
import matplotlib.pyplot as plt
import read
from scipy.io import FortranFile

###############################
#   read binary files
###############################
ixjx = read.read_grid2d('ixjx.dat')
ix = ixjx[0]
jx = ixjx[1]
margin = ixjx[2]
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
ay= read.read_phys2d('ay.dat',ix,jx)


##############
# setting
##############
nd = np.size(t)
x_1d = x[0,margin:]
z_1d = z[margin:,0]

x_ori_min = x[0,0]
xmin = x_1d[0]
dx = (xmin-x_ori_min)/margin
x_interval = np.arange(-xmin+dx,xmin,dx)
x_1d_all = np.concatenate([-x_1d[::-1],x_interval,x_1d])
x_all_len = np.size(x_1d_all)

z_ori_min = z[0,0]
zmin = z_1d[0]
dz = (zmin-z_ori_min)/margin
z_interval = np.arange(-zmin+dz,zmin,dz)
z_1d_all = np.concatenate([-z_1d[::-1],z_interval,z_1d])
z_all_len = np.size(z_1d_all)


## it is possible to solve z=0, so z_interval should be empty 
if(np.size(z_interval) != 0):
   exit() 


X, Z = np.meshgrid(x_1d_all,z_1d_all)

ro_all = np.empty((x_all_len,z_all_len))
ro_all[:,:] = np.nan













## for draw graph
fig = plt.figure(figsize=(8,8))
plt.rcParams['font.size']=10
plt.rcParams['font.family']='STIXGeneral'
plt.rcParams['mathtext.fontset']='stix'

ax = fig.add_subplot(111)
ax.set_aspect('equal')



n = 0
for n in range(nd):
    ax.set_xlim(0,1.5)
    ax.set_ylim(0,1.5)
    tle = r'$t=$ '+str(t[n])[:7]+'/ Density, Poloidal Magnetic Field line'
    #im1 = ax.pcolormesh(x,z,np.log10(ro[n]),cmap='jet',vmin=0,vmax=3)
    im1 = ax.contourf(x, z, np.log10(ro[n]),
                       levels=np.linspace(-4,1,11),extend='both',
                       vmin=-4,vmax=1,cmap='jet')
    ax.set_title(r"Density $\log_{10}(\rho)$")
    #ax.contour(x,z,ay[n],levels=20,colors='white',linewidths=0.5)
    #ax.set_title(r"$\log_{10}(\rho)$")

    #im1 = ax.contourf(x, z, np.log10(pr[n]),
    #                   levels=np.linspace(-10,1,11),extend='both',
    #                   vmin=-4,vmax=1,cmap='jet')
    #ax.set_title(r"$\log_{10}(P_{\mathrm{gas}})$")
    #im1 = ax.pcolormesh(x,z,np.sqrt(vx[n]**2+vz[n]**2),cmap='plasma',vmin=0,vmax=0.3)
    #im1 = ax.pcolormesh(x,z,np.log(pr[n]),cmap='plasma',vmin=-8,vmax=0)
    #im1 = ax.pcolormesh(x,z,bz[n],cmap='plasma')
    im2 = ax.contour(x,z,ay[n],colors='w',levels=15)
    im2.clabel(fmt='%5.4f')
    #im2 = ax.streamplot(X,Z,bx[n],bz[n],density=0.75,color='w')
    cbar1 = fig.colorbar(im1,ax=ax)
    #ax.set_title(tle)

    #plt.show()
    #exit()

    savename = 'fig/'+str(n).zfill(3)+'.png'
    plt.savefig(savename,dpi=300)
    cbar1.remove()
    ax.clear()


exit()
