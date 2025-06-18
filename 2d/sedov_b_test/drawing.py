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
    im1 = ax.pcolormesh(x,z,vx[n],cmap='plasma')
    cbar1 = fig.colorbar(im1,ax=ax)

    savename = 'fig/'+str(n).zfill(3)+'.png'
    plt.savefig(savename,dpi=300)
    cbar1.remove()
    ax.clear()



exit()
ax1 = fig.add_subplot(221)
ax1.set_aspect('equal')
ax2 = fig.add_subplot(222)
ax2.set_aspect('equal')
ax3 = fig.add_subplot(223)
ax3.set_aspect('equal')

for n in range(0,nd):
    ax1.set_ylim(0,1)
    ax1.set_xlim(0,1)
    ax2.set_ylim(0,1)
    ax2.set_xlim(0,1)
    ax3.set_ylim(0,1)
    ax3.set_xlim(0,1)
    tle = str(t[n])[:6]
    ax1tle = ax1.set_title(r"$\rho$ (Density) : $t=$"+tle)
    ax2tle = ax2.set_title(r"$\log P_r$ (Pressure) : $t=$"+tle)
    ax3tle = ax3.set_title(r"$V_{\text{all}}$ (Velocity) : $t=$"+tle)
    im1 = ax1.pcolormesh(x,z,ro[n],vmin=0,vmax=3,cmap='plasma')
    im2 = ax2.pcolormesh(x,z,np.log(pr[n]),vmin=-8,vmax=-4,cmap='plasma')
    #im3 = ax3.pcolormesh(x,z,np.sqrt(vx[n]**2+vz[n]**2),vmin=0, vmax=0.1, cmap='plasma')
    cbar1 = fig.colorbar(im1,ax=ax1)
    cbar2 = fig.colorbar(im2,ax=ax2)
    #cbar3 = fig.colorbar(im3,ax=ax3)
    svnm = "fig/"+ str(n).zfill(3) + ".png"
    plt.savefig(svnm,dpi=300)
    cbar1.remove()
    cbar2.remove()
    #cbar3.remove()
    ax1.clear()
    ax2.clear()
    #ax3.clear()

#plt.show()

exit()


#for n in range(2):
#    for i in range(ix):
#        for j in range(jx):
#            print(ro[n][i][j], end=' ')
#        print()
#    print()
#
#exit()

print('x=')
print(x)
print('z=')
print(z)
print()

#print(x)
#print(z)
for n in range(5):
    print(pr[n])
    print()
#n = 1
#exit()

#print(t[:])
#for n in range(7,8):
#    for i in range(0,4):
#        for j in range(0,4):
#            print(pr[n][i][j],end=' ')
#        print()
#exit()


# for graphing
fig = plt.figure(figsize=(10, 5))
plt.rcParams['font.size']=12
plt.rcParams['font.family']='STIXGeneral'
plt.rcParams['mathtext.fontset']='stix'


#ax00 = fig.add_subplot(111,projection='3d')
##ax00.set_xlim(0,1.0)
##ax00.set_ylim(0,5)
###ax00.set_ylim(-0.2,1.2)
##ax00.set_zlim(0,1.2)
#ax00.set_zlim(1e-4,1)
#ax00.set_zscale('log')
#ax00.set_title(r'Density')
#n = 0
#ax00.plot_surface(x,z,pr[0])
#plt.show()

#exit()


ax00 = fig.add_subplot(221)
ax01 = fig.add_subplot(222)
ax01.set_xlim(0,1.0)
#ax01.set_ylim(-0.2,1.2)
ax01.set_ylim(1e-4,1)
ax01.set_yscale('log')
ax01.set_title(r'Pressure')

ax10 = fig.add_subplot(223)
ax10.set_xlim(0,1.0)
ax10.set_ylim(-0.05,0.15)
#ax10.set_ylim(-0.2,1.2)
ax10.set_title(r'Velocity ($x$)')

for i in range(0,np.size(t)-1,2):
    if (i==1):
        continue
    lb = r'$t=$'+str(t[i])[:4]

    ax00.plot(diag, diag_ro, label=lb)
    ax01.plot(diag, diag_pr, label=lb)
    ax10.plot(diag, diag_vx, label=lb)

ax00.legend(loc='upper left', bbox_to_anchor=(1,1))
ax01.legend(loc='upper left', bbox_to_anchor=(1,1))
ax10.legend(loc='upper left', bbox_to_anchor=(1,1))

plt.tight_layout()

#plt.savefig('sedov1d.png', dpi=300,bbox_inches='tight')
plt.show()

exit()

