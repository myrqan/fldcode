import numpy as np
import matplotlib.pyplot as plt
import read

t = read.read_0d('t.dac')
x = read.read_0d('x.dac')
vx = read.read_1d('vx.dac')
rho = read.read_1d('rho.dac')
p = read.read_1d('p.dac')
#print(max(vx[0]))


# for graphing
fig = plt.figure(figsize=(10, 5))
plt.rcParams['font.size']=12
plt.rcParams['font.family']='STIXGeneral'
plt.rcParams['mathtext.fontset']='stix'


ax00 = fig.add_subplot(221)
ax00.set_xlim(0,1.0)
ax00.set_ylim(0,10)
##ax00.set_ylim(-0.2,1.2)
ax00.set_title(r'Density')

ax01 = fig.add_subplot(222)
ax01.set_xlim(0,1.0)
#ax01.set_ylim(-0.2,1.2)
ax01.set_ylim(1e-4,1)
ax01.set_yscale('log')
ax01.set_title(r'Pressure')

ax10 = fig.add_subplot(223)
ax10.set_xlim(0,1.0)
ax10.set_ylim(-0.05,0.20)
#ax10.set_ylim(-0.2,1.2)
ax10.set_title(r'Velocity ($x$)')

for i in range(0,np.size(t)-1,2):
    if (i==1):
        continue
    lb = r'$t=$'+str(t[i])[:4]

    ax00.plot(x, rho[i],label=lb)
    ax01.plot(x, p[i],label=lb)
    ax10.plot(x, vx[i],label=lb)

ax00.legend(loc='upper left', bbox_to_anchor=(1,1))
ax01.legend(loc='upper left', bbox_to_anchor=(1,1))
ax10.legend(loc='upper left', bbox_to_anchor=(1,1))

plt.tight_layout()

plt.savefig('sedov1d.png', dpi=300,bbox_inches='tight')
plt.show()

exit()

