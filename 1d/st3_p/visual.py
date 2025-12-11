import numpy as np
import matplotlib.pyplot as plt

tt = np.fromfile('dat/tt.dat',dtype=np.float64)
xx = np.fromfile('dat/x.dat',dtype=np.float64)

ns = np.size(tt)
ix = np.size(xx)

ro = np.fromfile('dat/ro.dat',dtype=np.float64).reshape(ns,ix)
vx = np.fromfile('dat/vx.dat',dtype=np.float64).reshape(ns,ix)
vy = np.fromfile('dat/vy.dat',dtype=np.float64).reshape(ns,ix)
pr = np.fromfile('dat/pr.dat',dtype=np.float64).reshape(ns,ix)

n = ns-1
for i in range(ix):
    print(ro[n,i])
exit()







exit()

## for draw graph
fig = plt.figure(figsize=(7, 5))
plt.rcParams['font.size']=13
plt.rcParams['font.family']='STIXGeneral'
plt.rcParams['mathtext.fontset']='stix'
ax = fig.add_subplot(111)

#ax.set_xlim(0.45,0.55)



#for i in range(ix):
#    eps = 0.01
#    if(np.abs(xx[i]-0.75)<eps):
#        print(xx[i],ro[ns-2,i])
#
#
#exit()
#for i in range(0,ns,3):
n = ns-1






ax.plot(xx,ro[i])
ax.plot(xx,vx[i])
ax.plot(xx,pr[i])

plt.savefig("save0.png",dpi=360)
plt.show()




