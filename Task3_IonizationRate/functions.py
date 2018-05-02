import numpy as np

__all__ = ['find_nearest','fit_exp_linear']


def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]
    
def fit_exp_linear(x, y):
    y_log = np.log(y)
    K, A_log = np.polyfit(x, y_log,1)
    A = np.exp(A_log)
    return A, K
