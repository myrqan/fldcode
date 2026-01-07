import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os


### read parameter file
pname = 'param.txt'
targets=['mgn','gix','giz','mpx','mpz']
mgn = gix = giz = mpx = mpz = 0

with open(pname,'r',encoding='utf-8') as f:
    for line in f:
        if ':' in line:
            key,val = line.split(':',1)
            key = key.strip()
            if key in targets:
                val_str = val.strip()
                exec(f"{key}=int(float({val_str}))")


gx = gix-2*mgn
gz = giz-2*mgn
lx = gx//mpx
lz = gz//mpz


## load grid data
xx = np.fromfile('dat/x.dat',dtype=np.float64)
zz = np.fromfile('dat/z.dat',dtype=np.float64)


## load time data
try:
    times = np.fromfile('dat/t.dat',dtype=np.float64)
except FileNotFoundError:
    print('Wraning! dat/t.dat not found.')
    exit()


## load physical value
def load_data(var,steps,datdir='dat'):
    num_steps=len(steps)
    data3d = np.zeros((num_steps,gx,gz),dtype=np.float64)
    stimes = []

    for t_idx, step in enumerate(steps):
        if times is not None and step < len(times):
            t_val = times[step]
        else:
            t_val = float(step)
        stimes.append(t_val)

        for i in range(mpx):
            for j in range(mpz):
                fname= f"{datdir}/{step:03d}_{var}_{i}_{j}.dat"
                if not os.path.exists(fname):
                    continue
                raw = np.fromfile(fname,dtype=np.float64)
                ldata = raw.reshape((lx,lz),order='F')

                istart=i*lx
                iend = istart+lx
                jstart=j*lz
                jend = jstart+lz

                data3d[t_idx,istart:iend,jstart:jend] = ldata

    return data3d,np.array(stimes)



### changeable
target_steps = [0]
target_var = 'ro'

ro,time_vals = load_data('ro',target_steps)
#pr,time_vals = load_data('pr',target_steps)
#bz,time_vals = load_data('bz',target_steps)

x,z=np.meshgrid(xx,zz,indexing='ij')

## for draw graph
fig = plt.figure(figsize=(7, 5))
plt.rcParams['font.size']=11
plt.rcParams['font.family']='Avenir'
plt.rcParams['mathtext.fontset']='stix'
ax = fig.add_subplot(111)
ax.set_aspect('equal')

ax.set_xlim(xx.min(),xx.max())
#ax.set_ylim(0,3.0)
ax.set_ylim(zz.min(),zz.max())

## select frame index < len(targt_steps)
frame_idx = 0

current_time = time_vals[frame_idx]
current_step = target_steps[frame_idx]

cmap_obj=sns.color_palette("coolwarm",as_cmap=True)

cont = ax.contourf(
        x,
        z,
        ro[frame_idx,:,:],
        #np.log10(ro[frame_idx,:,:].T),
        #levels=np.linspace(-4,1,11),
        #extend='both',
        #vmin=-4,
        #vmax=1,
        cmap=cmap_obj
        )

cbar = plt.colorbar(cont,ax=ax)
cbar.set_label(r'$\rho$')

plt.show()
