import numpy as np
from scipy.io import FortranFile

def rd_reals0d(f):
    list = []
    while True:
        try:
            list.extend(f.read_reals(dtype='f8'))
        except:
            break
    return np.array(list)

def rd_reals1d(f):
    list = []
    while True:
        try:
            list.append(f.read_reals(dtype='f8'))
        except:
            break
    return np.array(list)



def read_0d(file):
    ## t, x
    f = FortranFile(file, 'r')
    list_0d = rd_reals0d(f)
    return list_0d


def read_1d(file):
    ## ro, pr, vx
    f = FortranFile(file, 'r')
    list_1d = rd_reals1d(f)
    return list_1d



