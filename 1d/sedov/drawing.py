import numpy as np
import matplotlib.pyplot as plt
import read

t = read.read_0d('t.dac')
x = read.read_0d('x.dac')
vx = read.read_1d('vx.dac')
rho = read.read_1d('rho.dac')
p = read.read_1d('p.dac')
#print(max(vx[0]))


def plotall():
# for graphing
    fig = plt.figure(figsize=(8.5, 5))
    plt.rcParams['font.size']=15
    plt.rcParams['font.family']='STIXGeneral'
    plt.rcParams['mathtext.fontset']='stix'
    plt.ylim(0,1)
    i = 0
    plt.plot(x, vx[i],label=r'$v_r$ (Velocity)')
    plt.plot(x, p[i],label=r'$p$ (Pressure)')
    plt.plot(x,rho[i],label=r'$\rho$ (Density)')
    graph_title = r'$t = $' + str(t[i])[:4]
    plt.title(graph_title)
    plt.legend(loc='upper left', bbox_to_anchor=(1,1))
    plt.tight_layout()
    plt.show()
    plt.cla()

def plotr():
    # for graphing
    fig = plt.figure(figsize=(8.5, 5))
    plt.rcParams['font.size']=15
    plt.rcParams['font.family']='Palatino'
    plt.rcParams['mathtext.fontset']='stix'
    ax = fig.add_subplot(111)
    plt.xlim(0,1)
    for i in range(1, np.size(t)):
        lb = r'$t=$'+str(t[i])[:4]
        ax.plot(x,rho[i],label=lb)
    plt.legend(loc='upper left', bbox_to_anchor=(1,1))
    plt.tight_layout()
    plt.show()
    plt.cla()


# for graphing
fig = plt.figure(figsize=(10, 5))
plt.rcParams['font.size']=12
plt.rcParams['font.family']='STIXGeneral'
plt.rcParams['mathtext.fontset']='stix'


ax00 = fig.add_subplot(221)
ax00.set_xlim(0,1.0)
ax00.set_ylim(0,5)
ax00.set_title(r'Density')

ax01 = fig.add_subplot(222)
ax01.set_xlim(0,1.0)
#ax01.set_ylim(0,3)
ax01.set_ylim(1e-4,1)
ax01.set_yscale('log')
ax01.set_title(r'Pressure')

ax10 = fig.add_subplot(223)
ax10.set_xlim(0,1.0)
ax10.set_ylim(-0.05,0.10)
#ax10.set_ylim(0,1)
ax10.set_title(r'Velocity ($x$)')

for i in range(0,np.size(t), 3):
    lb = r'$t=$'+str(t[i])[:4]

    ax00.plot(x, rho[i],label=lb)
    ax01.plot(x, p[i],label=lb)
    ax10.plot(x, vx[i],label=lb)

ax00.legend(loc='upper left', bbox_to_anchor=(1,1))
ax01.legend(loc='upper left', bbox_to_anchor=(1,1))
ax10.legend(loc='upper left', bbox_to_anchor=(1,1))

plt.tight_layout()
plt.show()
exit()

