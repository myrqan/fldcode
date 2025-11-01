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





## for draw graph
fig = plt.figure(figsize=(7, 5))
plt.rcParams['font.size']=13
plt.rcParams['font.family']='STIXGeneral'
plt.rcParams['mathtext.fontset']='stix'
ax = fig.add_subplot(111)

for n in [0,5,10,14]:
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





