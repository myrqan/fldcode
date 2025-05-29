import numpy as np
from scipy.io import FortranFile

def read_prop(f):
    byo = f.read_record('i4')
    ver = f.read_record('i4')
    typ = f.read_record('i4')
    dim = f.read_record('i4')
    grs = f.read_record('i4')
    if (dim == 3):
        return [grs[0],grs[1]]
    return(grs)

def read_t(file):
    f = FortranFile(file,'r')
    #read_prop(f)
    t = []
    while True:
        try:
            t.extend(f.read_record('f8'))
        except:
            break
    return np.array(t)

def read_x2d(file,ix,jx):
    f = FortranFile(file,'r')
    #read_prop(f)
    x = f.read_record('f8').reshape(ix,jx,order='F')
    return np.array(x)

def read_phys2d(file,ix,jx):
    f = FortranFile(file,'r')
    #list = read_prop(f)
    #ix = list[0]
    #jx = list[1]
    phys = []
    while True:
        try:
            phys.append(f.read_record('f8').reshape(ix,jx,order='F'))
        except:
            break
    return np.array(phys)


def read_grid2d(file):
    f = FortranFile(file,'r')
    ixjx = f.read_record('i4')
    return np.array(ixjx)

