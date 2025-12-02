import numpy as np
import matplotlib.pyplot as plt
import read
import seaborn as sns

datdir = 'dat/'

# --- データ読み込み ---
tt = read.rd_0d(datdir+'t.dat')
xx = read.rd_0d(datdir+'x.dat')
uu = read.rd_1d(datdir+'uu.dat')
# ro = read.rd_1d(datdir+'ro.dat')
# vx = read.rd_1d(datdir+'vx.dat')
# vy = read.rd_1d(datdir+'vy.dat')
# bx = read.rd_1d(datdir+'bx.dat')
# by = read.rd_1d(datdir+'by.dat')
# pr = read.rd_1d(datdir+'pr.dat')
# ei = read.rd_1d(datdir+'ei.dat')

## for draw graph
fig = plt.figure(figsize=(7, 9))
#fig, (ax,ax2) = plt.subplots(2,1,figsize=(8,6),sharex=True)
plt.rcParams['font.size']=13
plt.rcParams['font.family']='STIXGeneral'
plt.rcParams['mathtext.fontset']='stix'
ax = fig.add_subplot(211)
ax2 = fig.add_subplot(212,sharex=ax)
ax.set_xlim(0,2)
ax2.ticklabel_format(style='sci',axis='y',scilimits=(0,0))


ax.set_title(r"numerical ($+$) v.s. analytical solution ($-$)")
ax2.set_title(rf"numerical $-$ analytical; gx $=$ {np.size(xx)}")


num_steps = 5
#colors = plt.cm.plasma(np.linspace(0, 1, num_steps))
colors = sns.color_palette("deep",n_colors=num_steps)

m_sp = int(np.size(xx)/100)

w = 0.05
mid = 0.5
eta = 0.1
aa = 1.0
#aa = 0.0

qq = np.zeros_like(uu)
for i in range(np.size(tt)):
    qq[i] = w/np.sqrt(w**2+4*eta*tt[i])*np.exp(-(xx-(mid+aa*tt[i]))**2/(w**2+4*eta*tt[i]))

for i in range(num_steps):
    tle = r"$t=$"+ str(tt[i])[:5]
    ax.plot(xx,uu[i],'+',markevery=m_sp//2,label=tle,color=colors[i])
    ax.plot(xx,qq[i],linestyle='-',color=colors[i])
    ax2.plot(xx,(uu[i]-qq[i]),label=tle)
#print(xx)

plt.tight_layout()
plt.legend()
plt.savefig(f"gx{np.size(xx)}eta1e-1.png",dpi=300)
#plt.show()
