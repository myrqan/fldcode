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
xmin = x_1d[0]
zmin = z_1d[0]

zz_1d = np.append(-z_1d[::-1],z_1d)
xx, zz = np.meshgrid(x_1d, zz_1d)


ro = ro[:,margin:,margin:]
vx = vx[:,margin:,margin:]
vz = vz[:,margin:,margin:]
by = by[:,margin:,margin:]
pr = pr[:,margin:,margin:]
ay = ay[:,margin:,margin:]

roz = np.append(ro[:,::-1,:],ro,axis=1)
vxz = np.append(vx[:,::-1,:],vx,axis=1)
vzz = np.append(-vz[:,::-1,:],vz,axis=1)
byz = np.append(-by[:,::-1,:],by,axis=1)
prz = np.append(pr[:,::-1,:],pr,axis=1)
ayz = np.append(ay[:,::-1,:],ay,axis=1)


## for draw graph
fig = plt.figure(figsize=(5.5, 9))
plt.rcParams['font.size']=13
plt.rcParams['font.family']='STIXGeneral'
plt.rcParams['mathtext.fontset']='stix'
ax = fig.add_subplot(111)

## flag: density, bphi, pressure, beta
flag = 'pressure'


for n in range(nd):
    ax.set_aspect('equal')
    ax.set_xlim(0,0.5)
    ax.set_ylim(-0.5,0.5)
    tle = r'$t=$ '+str(t[n])[:5]
    if(flag == 'density'):
        im1 = ax.contourf(xx,zz, np.log10(roz[n]),
                           levels=np.linspace(-4,1,11),extend='both',
                           vmin=-4,vmax=1,cmap='jet')
        cbar1 = fig.colorbar(im1,ax=ax)
        ax.set_title(tle+r" Density $\log_{10}(\rho)$")
        ax.contour(xx, zz, np.log10(roz[n]),
                           levels=np.linspace(-4,1,31),extend='both',
                           vmin=-4,vmax=1,colors='black',
                           linewidths=0.8)
    elif(flag=='bphi'):
        im1 = ax.contourf(xx,zz,byz[n],
                           levels=np.linspace(-3,3,11),extend='both',
                           vmin=-3,vmax=3,cmap='jet')
        cbar1 = fig.colorbar(im1,ax=ax)
        ax.set_title(tle+r" toroidal magnetic field $B_{\phi}$")
        ax.contour(xx, zz, byz[n],
                           levels=np.linspace(-3,3,31),extend='both',
                           vmin=-3,vmax=3,colors='black',
                           linewidths=0.8)
    elif(flag=='pressure'):
        im1 = ax.contourf(xx,zz,np.log10(prz[n]),
                           levels=np.linspace(-3,2,11),extend='both',
                           vmin=-3,vmax=2,cmap='jet')
        cbar1 = fig.colorbar(im1,ax=ax)
        ax.set_title(tle+r" gas pressure $\log_{10}(P_{\mathrm{gas}})$")
        ax.contour(xx, zz, np.log10(prz[n]),
                           levels=np.linspace(-3,2,31),extend='both',
                           vmin=-3,vmax=2,colors='black',
                           linewidths=0.8)

    ax.quiver(xx[::5,::5], zz[::5,::5], vxz[n,::5,::5], vzz[n,::5,::5],
              scale=10,color='black')

    #im1 = ax.contourf(x, z, np.log10(pr[n]),
    #                   levels=np.linspace(-10,1,11),extend='both',
    #                   vmin=-4,vmax=1,cmap='jet')
    #ax.set_title(r"$\log_{10}(P_{\mathrm{gas}})$")
    #im1 = ax.pcolormesh(x,z,np.sqrt(vx[n]**2+vz[n]**2),cmap='plasma',vmin=0,vmax=0.3)
    #im1 = ax.pcolormesh(x,z,np.log(pr[n]),cmap='plasma',vmin=-8,vmax=0)
    #im1 = ax.pcolormesh(x,z,bz[n],cmap='plasma')
    im2 = ax.contour(xx,zz,ayz[n],colors='w',levels=np.arange(0.01,0.5,0.01))
    #im2.clabel(fmt='%5.4f')
    #ax.set_title(tle)

    #plt.show()
    #exit()

    savename = 'fig/'+str(n).zfill(3)+'.png'
    plt.savefig(savename,dpi=300)
    cbar1.remove()
    ax.clear()

exit()

