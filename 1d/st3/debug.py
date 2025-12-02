import numpy as np
import matplotlib.pyplot as plt
import read

datdir = 'dat/'

tt = read.rd_0d(datdir+'t.dat')
xx = read.rd_0d(datdir+'x.dat')
ro = read.rd_1d(datdir+'ro.dat')
#vx = read.rd_1d(datdir+'vx.dat')
#pr = read.rd_1d(datdir+'pr.dat')
#ei = read.rd_1d(datdir+'ei.dat')

# be consistent
ix = np.size(xx)
xx = np.delete(xx,[0,ix-1],0)
ro = np.delete(ro,[0,ix-1],1)


ns = np.size(tt)
ix = np.size(xx)

print(ns,ix)

for n in range(0,10):
    for i in range(ix):
        if(xx[i] >= 0.45 and xx[i] <= 0.55):
            print(xx[i],ro[n,i])
    print('---')

