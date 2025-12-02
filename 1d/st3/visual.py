import numpy as np
import matplotlib.pyplot as plt
import read

datdir = 'dat/'

tt = read.rd_0d(datdir+'t.dat')
xx = read.rd_0d(datdir+'x.dat')
ro = read.rd_1d(datdir+'ro.dat')
vx = read.rd_1d(datdir+'vx.dat')
pr = read.rd_1d(datdir+'pr.dat')
ei = read.rd_1d(datdir+'ei.dat')


ns = np.size(tt)


## for draw graph
fig = plt.figure(figsize=(7, 5),)
plt.rcParams['font.size']=13
plt.rcParams['font.family']='STIXGeneral'
plt.rcParams['mathtext.fontset']='stix'
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)


#ax.plot(xx,ro[10])
#plt.show()


ax1.contourf(xx,tt,ro)
ax2.contourf(xx,tt,vx)
ax3.contourf(xx,tt,pr)
for n in range(0,ns,ns//3):
    ax4.plot(xx,ro[n])
plt.show()


exit()

for n in [0,5,10]:
    ax.set_xlim(0,1)
    time = tt[n]
    ax.set_title(r'$t =$ '+str(time)[:5])

    ax.plot(xx,ro[n],label=r'density')
    ax.plot(xx,vx[n],label=r'velocity')
    ax.plot(xx,ei[n],label=r'internal energy')

    plt.legend()
    savname = 'graph/' + str(n).zfill(3) + '.png'
    plt.savefig(savname,dpi=300)
    plt.cla()


#plt.show()





