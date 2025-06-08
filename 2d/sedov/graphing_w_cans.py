import numpy as np
import matplotlib.pyplot as plt
import read
from scipy.io import FortranFile
import seaborn as sns

def read_prop_cans(f):
    byo = f.read_record('i4')
    ver = f.read_record('i4')
    typ = f.read_record('i4')
    dim = f.read_record('i4')
    grs = f.read_record('i4')
    if (dim == 3):
        return [grs[0],grs[1]]
    return(grs)

def read_t_cans(file):
    f = FortranFile(file,'r')
    read_prop_cans(f)
    t = []
    while True:
        try:
            t.extend(f.read_record('f8'))
        except:
            break
    return np.array(t)

def read_x_cans(file):
    f = FortranFile(file,'r')
    read_prop_cans(f)
    x = np.array(f.read_record('f8'))
    return x

def read_phys2d_cans(file):
    f = FortranFile(file,'r')
    list = read_prop_cans(f)
    ix = list[0]
    jx = list[1]
    phys = []
    while True:
        try:
            phys.append(f.read_record('f8').reshape(ix,jx,order='F'))
        except:
            break
    return np.array(phys)

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

#### read cans output
dir = '/Users/crutont/Documents/cans/cans2d/md_sedov/hdmlw/'
t_cans = read_t_cans(dir+'t.dac')
x_cans = read_x_cans(dir+'x.dac')
z_cans = read_x_cans(dir+'z.dac')

ix_cans = np.size(x_cans)
jx_cans = np.size(z_cans)

X_cans,Z_cans = np.meshgrid(x_cans,z_cans)

ro_cans = read_phys2d_cans(dir+'ro.dac')
pr_cans = read_phys2d_cans(dir+'pr.dac')
vx_cans = read_phys2d_cans(dir+'vx.dac')
vz_cans = read_phys2d_cans(dir+'vz.dac')


if(ix != ix_cans):
    print(ix, ix_cans)
    print('Error!: grid size is different')
    exit()

nd = np.size(t)
nd_cans = np.size(t_cans)


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

diag_cans = np.zeros(ix)
diag_ro_cans = np.zeros((nd_cans,ix))
diag_pr_cans = np.zeros((nd_cans,ix))
diag_vx_cans = np.zeros((nd_cans,ix))
diag_vz_cans = np.zeros((nd_cans,ix))
for i in range(ix):
    diag_cans[i] = x_cans[i]

for n in range(nd):
    diag_ro_cans[n] = np.diag(ro_cans[n])
    diag_pr_cans[n] = np.diag(pr_cans[n])
    diag_vx_cans[n] = np.diag(vx_cans[n])
    diag_vz_cans[n] = np.diag(vz_cans[n])


# for graphing
fig = plt.figure(figsize=(15, 8))
plt.rcParams['font.size']=12
plt.rcParams['font.family']='STIXGeneral'
plt.rcParams['mathtext.fontset']='stix'
colors = sns.color_palette("muted", nd).as_hex()

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
    ax00.plot(diag, diag_ro[i], label=lb,color=colors[i%len(colors)])
    ax00.plot(diag_cans, diag_ro_cans[i],'--',color=colors[i%len(colors)])

    ax01.plot(diag, diag_pr[i], label=lb,color=colors[i%len(colors)])
    ax01.plot(diag_cans, diag_pr_cans[i],'--',color=colors[i%len(colors)])

    ax10.plot(diag, np.sqrt(diag_vx[i]**2+diag_vz[i]**2), label=lb,color=colors[i%len(colors)])
    ax10.plot(diag_cans, np.sqrt(diag_vx_cans[i]**2+diag_vz_cans[i]**2), '--',color=colors[i%len(colors)])

ax00.legend(loc='upper left', bbox_to_anchor=(1,1))
ax01.legend(loc='upper left', bbox_to_anchor=(1,1))
ax10.legend(loc='upper left', bbox_to_anchor=(1,1))

plt.tight_layout()

plt.savefig('fig/sedov2d_w_cans.png',dpi=300,bbox_inches='tight')
#plt.show()

exit()


