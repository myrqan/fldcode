import numpy as np
import matplotlib.pyplot as plt
import read
from scipy.io import FortranFile
import seaborn as sns

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
vy = vy[:,margin:,margin:]
vz = vz[:,margin:,margin:]
bx = bx[:,margin:,margin:]
by = by[:,margin:,margin:]
bz = bz[:,margin:,margin:]
pr = pr[:,margin:,margin:]
ay = ay[:,margin:,margin:]

roz = np.append(ro[:,::-1,:],ro,axis=1)
vxz = np.append(vx[:,::-1,:],vx,axis=1)
vyz = np.append(vy[:,::-1,:],vy,axis=1)
vzz = np.append(-vz[:,::-1,:],vz,axis=1)
bxz = np.append(bx[:,::-1,:],bx,axis=1)
byz = np.append(-by[:,::-1,:],by,axis=1)
bzz = np.append(bz[:,::-1,:],bz,axis=1)
prz = np.append(pr[:,::-1,:],pr,axis=1)
ayz = np.append(ay[:,::-1,:],ay,axis=1)

prmgz = 0.5*(bxz**2 + byz**2 + bzz**2)


## for draw graph
fig = plt.figure(figsize=(5.5, 9))
#fig = plt.figure(figsize=(5.5, 7))

plt.rcParams['font.size']=13
plt.rcParams['font.family']='STIXGeneral'
plt.rcParams['mathtext.fontset']='stix'
ax = fig.add_subplot(111)

## flag: density, bphi, pressure, beta
flag = 'density'


cmap_obj = sns.color_palette("coolwarm",as_cmap=True)

for n in range(nd):
    #if(n%10!=0):
        #continue
    ax.set_aspect('equal')
    ax.set_xlim(0,1.5)
    ax.set_ylim(-0.5,3.0)
    tle = r'$t=$ '+str(t[n])[:5]
    if(flag == 'density'):
        im1 = ax.contourf(xx,zz, np.log10(roz[n]),
                           levels=np.linspace(-4,1,11),extend='both',
                           vmin=-4,vmax=1,cmap=cmap_obj)
        cbar1 = fig.colorbar(im1,ax=ax)
        ax.set_title(tle+r" Density $\log_{10}(\rho)$")
        #ax.contour(xx, zz, np.log10(roz[n]),
        #                   levels=np.linspace(-4,1,31),extend='both',
        #                   vmin=-4,vmax=1,colors='black',
        #                   linewidths=0.5)
    elif(flag=='vy'):
        im1 = ax.contourf(xx,zz,vyz[n],
                           levels=np.linspace(-2,2,21),extend='both',
                           vmin=-2,vmax=2,cmap=cmap_obj)
        cbar1 = fig.colorbar(im1,ax=ax)
        ax.set_title(tle+r" velocity field $V_{\phi}$")
        #ax.contour(xx, zz, byz[n],
        #                  levels=np.linspace(-3,3,31),extend='both',
        #                  vmin=-3,vmax=3,colors='black',
        #                  linewidths=0.5)
    elif(flag=='bphi'):
        im1 = ax.contourf(xx,zz,byz[n],
                           levels=np.linspace(-0.5/np.sqrt(4*np.pi),0.5/np.sqrt(4*np.pi),21),extend='both',
                           vmin=-0.5/np.sqrt(4*np.pi),vmax=0.5/np.sqrt(4*np.pi),cmap=cmap_obj)
        cbar1 = fig.colorbar(im1,ax=ax)
        ax.set_title(tle+r" toroidal magnetic field $B_{\phi}$")
        #ax.contour(xx, zz, byz[n],
        #                   levels=np.linspace(-0.5,0.5,31),extend='both',
        #                   vmin=-0.5,vmax=0.5,colors='black',
        #                   linewidths=0.5)
    elif(flag=='pressure'):
        im1 = ax.contourf(xx,zz,np.log10(prz[n]),
                           levels=np.linspace(-3,2,11),extend='both',
                           vmin=-3,vmax=2,cmap=cmap_obj)
        cbar1 = fig.colorbar(im1,ax=ax)
        ax.set_title(tle+r" gas pressure $\log_{10}(P_{\mathrm{gas}})$")
        ax.contour(xx, zz, np.log10(prz[n]),
                           levels=np.linspace(-3,2,31),extend='both',
                           vmin=-3,vmax=2,colors='black',
                           linewidths=0.5)

    elif(flag=='beta'):
        im1 = ax.contourf(xx,zz,np.log10(prz[n]/prmgz[n]),
                           levels=np.linspace(-3,3,17),extend='both',
                           vmin=-3,vmax=3,cmap=cmap_obj)
        cbar1 = fig.colorbar(im1,ax=ax)
        ax.set_title(tle+r" plasma beta $\log_{10}(\beta)$")
        #ax.contour(xx, zz, np.log10(prz[n]/prmgz[z]),
                           #levels=np.linspace(-3,2,31),extend='both',
                           #vmin=-3,vmax=2,colors='black',
                           #linewidths=0.5)

    ax.quiver(xx[::30,::30], zz[::30,::30], vxz[n,::30,::30], vzz[n,::30,::30],
              scale=1.,color='black')

    #im1 = ax.contourf(x, z, np.log10(pr[n]),
    #                   levels=np.linspace(-10,1,11),extend='both',
    #                   vmin=-4,vmax=1,cmap='jet')
    #ax.set_title(r"$\log_{10}(P_{\mathrm{gas}})$")
    #im1 = ax.pcolormesh(x,z,np.sqrt(vx[n]**2+vz[n]**2),cmap='plasma',vmin=0,vmax=0.3)
    #im1 = ax.pcolormesh(x,z,np.log(pr[n]),cmap='plasma',vmin=-8,vmax=0)
    #im1 = ax.pcolormesh(x,z,bz[n],cmap='plasma')
    im2 = ax.contour(xx,zz,ayz[n],colors='w',levels=np.arange(0.0,0.5,0.002),linewidths=0.6)
    
    #ax.set_title(tle)
    #im2 = ax.contour(xx,zz,ayz[n],colors='k',levels=np.arange(0,0.5,0.001),linewidths=0.9)
    #im2.clabel(fmt='%5.4f')
    #ax.set_title(tle)

    #plt.show()
    #exit()

    #savename = 'fig/'+str(n).zfill(3)+'.png'
    savename = 'fig/'+str(n).zfill(3)+'.png'
    plt.savefig(savename,dpi=300)
    cbar1.remove()
    ax.clear()
    if(n%50==0):
        print('png No. '+str(n)+' is completed.')

exit()

