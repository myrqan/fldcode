import numpy as np
import matplotlib.pyplot as plt
import read
from scipy.io import FortranFile

#t = read.read_0d('t.dac')
#x = read.read_0d('x.dac')
#vx = read.read_1d('vx.dac')
#rho = read.read_1d('rho.dac')
#p = read.read_1d('p.dac')
#print(max(vx[0]))


ixjx = read.read_grid2d('ixjx.dac')
ix = ixjx[0]
jx = ixjx[1]
t = read.read_t('t.dac')
x = read.read_x2d('x.dac',ix,jx)
z = read.read_x2d('z.dac',ix,jx)
ro= read.read_phys2d('rho.dac',ix,jx)
vx= read.read_phys2d('vx.dac',ix,jx)
vz= read.read_phys2d('vz.dac',ix,jx)
pr= read.read_phys2d('p.dac',ix,jx)

nd = np.size(t)

diag = np.zeros(ix)
diag_ro = np.zeros((nd,ix))
diag_pr = np.zeros((nd,ix))
diag_vx = np.zeros((nd,ix))
diag_vz = np.zeros((nd,ix))

for i in range(ix):
    diag[i] = x[i][i]

for n in range(nd):
    diag_ro[n] = np.diag(ro[n])
    diag_pr[n] = np.diag(pr[n])
    diag_vx[n] = np.diag(vx[n])
    diag_vz[n] = np.diag(vz[n])

# for graphing
fig = plt.figure(figsize=(15, 8))
plt.rcParams['font.size']=12
plt.rcParams['font.family']='STIXGeneral'
plt.rcParams['mathtext.fontset']='stix'


ax00 = fig.add_subplot(221)
ax01 = fig.add_subplot(222)
ax10 = fig.add_subplot(223)

ax00.set_xlim(0,1.0)
ax00.set_ylim(0,5)
ax00.set_title(r'Density')


ax01.set_xlim(0,1.0)
ax01.set_ylim(1e-4,1)
ax01.set_yscale('log')
ax01.set_title(r'Pressure')

ax10.set_xlim(0,1.0)
ax10.set_ylim(-0.05,0.15)
ax10.set_title(r'Velocity (all)')

for i in range(0,nd,2):
    lb = r'$t=$'+str(t[i])[:4]

    ax00.plot(diag, diag_ro[i], label=lb)
    ax01.plot(diag, diag_pr[i], label=lb)
    ax10.plot(diag, np.sqrt(diag_vx[i]**2+diag_vz[i]**2), label=lb)

ax00.legend(loc='upper left', bbox_to_anchor=(1,1))
ax01.legend(loc='upper left', bbox_to_anchor=(1,1))
ax10.legend(loc='upper left', bbox_to_anchor=(1,1))

plt.tight_layout()

plt.savefig('fig/sedov2d.png',dpi=300,bbox_inches='tight')
#plt.show()

exit()

