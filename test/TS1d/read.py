import numpy as np
from scipy.io import FortranFile

def rd_r0(f):
    list = []
    while True:
        try:
            list.extend(f.read_reals(dtype='f8'))
        except:
            break
    return np.array(list)


def rd_r1(f):
    list = []
    while True:
        try:
            list.append(f.read_reals(dtype='f8'))
        except:
            break
    return np.array(list)

def rd_0d(fname):
    f = FortranFile(fname,'r')
    qq = rd_r0(f)
    return qq

def rd_1d(fname):
    f = FortranFile(fname,'r')
    qq = rd_r1(f)
    return qq

