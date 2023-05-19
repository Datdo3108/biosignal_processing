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

def convolution(X, W):
    X = np.reshape(X, (-1,1))
    W = np.reshape(W, (-1,1))

    C = np.concatenate([np.zeros(len(W) - 1).reshape(-1,1), X, np.zeros(len(W) - 1).reshape(-1,1)]).reshape((1,-1))
    C_remove = np.intc(len(X) + len(W) - 1)

    C_row = np.concatenate([W, np.zeros(C_remove).reshape(-1,1)])
    C_mat = np.tile(C_row, (C_remove, 1))
    C_mat = C_mat.ravel()[:-C_remove]
    C_mat = C_mat.reshape((C_remove, -1))

    Y = np.matmul(C_mat, C.T)/np.sum(W)
    Y_start = np.intc(len(W)/2 - 1)
    Y_stop = np.intc(len(X) + Y_start)

    return Y[Y_start:Y_stop]