import numpy as np
import matplotlib.pyplot as plt

### read parameter file 
filename = 'param.txt'
targets = ['mgn', 'gix', 'giz', 'mpx', 'mpz']
with open(filename, 'r', encoding='utf-8') as f:
    for line in f:
        if ':' in line:
            key, val = line.split(':', 1)
            key = key.strip()
            if key in targets:
                val_str = val.strip()
                exec(f"{key} = int(float({val_str}))")


# gx = gix-2*mgn
# gz = giz-2*mgn
# lx = gx//mpx
# lz = gz//mpz

### read data

xx = np.fromfile('dat/x.dat',dtype=np.float64)
zz = np.fromfile('dat/z.dat',dtype=np.float64)


print(xx)
print(zz)




exit()

x,z = np.meshgrid(xx,zz)

phi = np.exp(-(x**2+z**2))




## for draw graph
fig = plt.figure(figsize=(7, 5))
plt.rcParams['font.size']=11
plt.rcParams['font.family']='Avenir'
plt.rcParams['mathtext.fontset']='stixsans'
ax = fig.add_subplot(111)


ax.contourf(x,z,phi)

plt.show()


