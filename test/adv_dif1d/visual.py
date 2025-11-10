import numpy as np
import matplotlib.pyplot as plt
import read

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
fig = plt.figure(figsize=(7, 5))
plt.rcParams['font.size']=13
plt.rcParams['font.family']='STIXGeneral'
plt.rcParams['mathtext.fontset']='stix'
ax = fig.add_subplot(111)



for i in range(3):
    tle = r"$t=$"+ str(tt[i])[:4]
    ax.plot(xx,uu[i],label=tle)
#print(xx)

plt.legend()
plt.show()
