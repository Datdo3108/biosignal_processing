import numpy as np

def get_fourier_coefficient(X, Y, N_freq, T):
    A = []
    B = []
    for i in range(len(N_freq)):
        if i == 0:
            a = np.mean(Y)
            b = 0
        else:
            a = np.mean(Y*np.cos(N_freq[i]*np.pi*X/T)*2)
            b = np.mean(Y*np.sin(N_freq[i]*np.pi*X/T)*2)

        A.append(a)
        B.append(b)
    return A, B

def fourier_series_spectrum(X, Y, N_freq, T):
    A, B = get_fourier_coefficient(X=X, Y=Y, N_freq=N_freq, T=T)
    F = []

    for i in range(len(N_freq)):
        f = np.sqrt(pow(A[i],2) + pow(B[i],2))
        F.append(f)
    return F

def fourier_series(X, Y, N_freq, T):
    A, B = get_fourier_coefficient(X=X, Y=Y, N_freq=N_freq, T=T)
    F = []

    for x in X:
        f = 0
        for i in range(len(N_freq)):
            f = f + A[i]*np.cos(i*np.pi*x/T) + B[i]*np.sin(i*np.pi*x/T)
        F.append(f)
    return F
