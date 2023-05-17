import numpy as np

def windowing(X,W):
    Y = np.multiply(X,W)

    return Y

def rect(n):
    W =  np.ones(n)*0.8
    return W

def hann(n):
    N = np.arange(n)
    W = np.power(np.sin(np.pi*N/(n-1)),2)

    return W

def blackman(n, alpha):
    N = np.arange(n)
    W = (1 - alpha)/2 - 1/2*np.cos(2*np.pi*N/(n-1)) + alpha/2*np.cos(4*np.pi*N/(n-1))

    return W

def hamming(n):
    a = 25/46
    N = np.arange(n)
    W = a - (1 - a)*np.cos(2*np.pi*N/(n-1))

    return W