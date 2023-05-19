import numpy as np

def get_fourier_coefficient(X, Y, N_freq, T):
    X = np.reshape(X, (-1,1))
    Y = np.reshape(Y, (-1,1))
    N_freq = np.reshape(N_freq, (-1,1))

    A = np.mean(np.cos(N_freq*np.pi*np.transpose(X)/T)*2*np.transpose(Y), axis=1)
    B = np.mean(np.sin(N_freq*np.pi*np.transpose(X)/T)*2*np.transpose(Y), axis=1)

    A[0] = np.mean(Y)
    B[0] = 0

    return A, B

def fourier_series_spectrum(X, Y, N_freq, T):
    A, B = get_fourier_coefficient(X=X, Y=Y, N_freq=N_freq, T=T)
    F = np.sqrt(np.power(A,2) + np.power(B,2))
    return F

def fourier_series(X, Y, N_freq, T):
    X = np.reshape(X, (-1,1))
    Y = np.reshape(Y, (-1,1))
    N_freq = np.reshape(N_freq, (-1,1))

    A, B = get_fourier_coefficient(X=X, Y=Y, N_freq=N_freq, T=T)

    F = np.sum(np.cos(X*np.transpose(N_freq)*np.pi/T)*A + np.sin(X*np.transpose(N_freq)*np.pi/T)*B, axis=1)

    return F

def time_frequency_domain(X, Y, N_freq, T):
    X = np.reshape(X, (-1,1))
    Y = np.reshape(Y, (-1,1))
    N_freq = np.reshape(N_freq, (-1,1))

    A, B = get_fourier_coefficient(X=X, Y=Y, N_freq=N_freq, T=T)

    A = np.reshape(A, (-1,1))
    B = np.reshape(B, (-1,1))
    
    TF = np.cos(N_freq*np.pi*np.transpose(X)/T)*A + np.sin(N_freq*np.pi*np.transpose(X)/T)*B
    return TF

def power_spectrum(X, Y, N_freq, T):
    A, B = get_fourier_coefficient(X=X, Y=Y, N_freq=N_freq, T=T)
    P = np.power(A,2) + np.power(B,2)

    return P
