import numpy as np
import matplotlib.pyplot as plt
from decimal import *
import read

x = read.read_0d('x.dac')
t = read.read_0d('t.dac')
ro = read.read_1d('rho.dac')
pr = read.read_1d('p.dac')
vx = read.read_1d('vx.dac')


X = np.size(x)
N = np.size(t)

gam = 5.0/3.0
#xi0 = 1.08
xi0=1.15

ro0 = 1.0
we = 0.1
enttl = (np.sqrt(np.pi)*we)**3/(gam-1)
scl = np.zeros(N)

for i in range (1, N):
    #scl[i] = xi0 * enttl**0.2 * t[i]**(0.35)
    scl[i] = xi0 * enttl**0.2 * t[i]**(0.40)

d_s = np.zeros(N)
d_s[:] = np.nan

for i in range(1, N):
    d_s[i] = 0.4 * scl[i] / t[i]

lam = np.zeros((N, X))
for i in range(1, N):
    for j in range(X):
        lam[i,j] = x[j] / scl[i]



sc_ro = np.zeros((N, X))
sc_pr = np.zeros((N, X))
sc_vx = np.zeros((N, X))

for i in range(1, N):
    for j in range(1, X):
        sc_pr[i, j] = pr[i,j]/2 *(gam+1) * d_s[i]**(-2)
        sc_ro[i, j] = ro[i,j]/(gam+1)*(gam-1)
        sc_vx[i, j] = vx[i,j]/2*(gam+1)/d_s[i]



fig = plt.figure(figsize=(10, 5))
plt.rcParams['font.size']=15
plt.rcParams['font.family']='STIXGeneral'
plt.rcParams['mathtext.fontset']='stix'

def plotr():
    ax00 = fig.add_subplot(221)
    ax00.set_xlim(0.0, 1.2)
    ax00.set_title(r'$\rho$, (density) scaled')
    for n in range(0,N,4):
        if n == 0:
            continue
        time = Decimal(str(t[n])).quantize(Decimal('0.01'),ROUND_HALF_UP)
        ax00.plot(lam[n], sc_ro[n], label=(r'$t=$'+str(time)))
    ax00.legend(loc='upper left', bbox_to_anchor=(1,1))
    #plt.savefig('fig/ro_scaled.png',dpi=300,bbox_inches='tight')
    #plt.cla()
    #plt.show()


def plotv():
    ax01 = fig.add_subplot(222)
    plt.xlim(0.0, 1.2)
    plt.title(r'$v_r$, (velocity) scaled')
    plt.ylim(0,1.2)
    for n in range(0,N,4):
        if n == 0:
            continue
        time = Decimal(str(t[n])).quantize(Decimal('0.01'),ROUND_HALF_UP)
        ax01.plot(lam[n], sc_vx[n], label=(r'$t=$'+str(time)))
    ax01.legend(loc='upper left', bbox_to_anchor=(1,1))
    #plt.savefig('fig/vr_scaled.png',dpi=300, bbox_inches='tight')
    #plt.cla()
    #plt.show()

    
def plotp():
    ax02 = fig.add_subplot(223)
    plt.xlim(0.0, 1.2)
    plt.title(r'$P_r$, (pressure) scaled')
    #plt.yscale('log')
    plt.ylim(0,1.2)
    for n in range(0, N, 4):
        if n == 0:
            continue
        time = Decimal(str(t[n])).quantize(Decimal('0.01'),ROUND_HALF_UP)
        ax02.plot(lam[n], sc_pr[n], label=(r'$t=$'+str(time)))
    ax02.legend(loc='upper left', bbox_to_anchor=(1,1))
    #plt.savefig('fig/pr_scaled.png',dpi=300,bbox_inches='tight')
    #plt.cla()
    #plt.show()


def plotrvp():
    ax03 = fig.add_subplot(224)
    plt.xlim(0.0, 1)
    #plt.title("all")
    #plt.title(r'$P_r$, (pressure) scaled')
    #plt.yscale('log')
    plt.ylim(0,1)
    for n in range(1, N, N-2):
        if n == 1:
            continue
        time = Decimal(str(t[n])).quantize(Decimal('0.01'),ROUND_HALF_UP)
        ax03.set_title(r"$t=$"+str(time))
        ax03.plot(lam[n], sc_pr[n], label="pressure")
        ax03.plot(lam[n], sc_ro[n], label="density")
        ax03.plot(lam[n], sc_vx[n], label="velocity")
    ax03.legend(loc='upper left', bbox_to_anchor=(1,1))
    #plt.savefig('fig/pr_scaled.png',dpi=300,bbox_inches='tight')
    #plt.cla()
    #plt.show()

plotv()
plotp()
plotr()
plotrvp()



plt.tight_layout()
plt.show()

#plt.savefig('fig/scaled.png',dpi=300,bbox_inches='tight')

exit()
