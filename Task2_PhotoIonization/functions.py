import numpy as np

__all__ = ['find_nearest']


def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]
