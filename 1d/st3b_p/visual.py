import numpy as np
import matplotlib.pyplot as plt

tt = np.fromfile('dat/tt.dat',dtype=np.float64)
xx = np.fromfile('dat/x.dat',dtype=np.float64)

ns = np.size(tt)
ix = np.size(xx)

ro = np.fromfile('dat/ro.dat',dtype=np.float64).reshape(ns,ix)
vx = np.fromfile('dat/vx.dat',dtype=np.float64).reshape(ns,ix)
vy = np.fromfile('dat/vy.dat',dtype=np.float64).reshape(ns,ix)
bx = np.fromfile('dat/bx.dat',dtype=np.float64).reshape(ns,ix)
by = np.fromfile('dat/by.dat',dtype=np.float64).reshape(ns,ix)
pr = np.fromfile('dat/pr.dat',dtype=np.float64).reshape(ns,ix)

## for draw graph

fig = plt.figure(figsize=(7, 5))
plt.rcParams['font.size']=13
plt.rcParams['font.family']='STIXGeneral'
plt.rcParams['mathtext.fontset']='stix'
ax = fig.add_subplot(111)
#ax.set_xlim(0.45,0.55)

i = ns-1
ax.plot(xx,ro[i],label=r'density')
ax.plot(xx,vx[i],label=r'$x$-velocity')
ax.plot(xx,vy[i],label=r'$y$-velocity')
ax.plot(xx,by[i]/np.sqrt(4*np.pi),label=r'$B_y/\sqrt{4\pi}$')
ax.plot(xx,pr[i],label='pressure')
plt.legend()
savname='graph/'+str(i).zfill(3)+'.png'





###### execute #####

plt.savefig(savname,dpi=360)
plt.close()



#
#    savname='graph/'+str(i).zfill(3)+'.png'
#
#    plt.savefig(savname,dpi=360)
#    plt.close()




