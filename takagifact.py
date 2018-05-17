import numpy as np

def takagi_for_unitary(A):
    ###################################################
    ### takes a unitary matrix A such that A^T = A, ###
    ###    returns unitary U such that A = U U^T    ###
    ###################################################

    N = A.shape[0]
    U = np.eye(N)

    while A.shape[0] > 0:   
        n = A.shape[0]

        # constructing a vector v such that A*ybar = y
        x = np.random.rand(n)+1j*np.random.rand(n)
        x /= np.sqrt(np.dot(x.conj(),x))
        Axbar = np.tensordot(A,x.conj(),(1,0))
        mu = np.dot(x.conj(),Axbar)

        if np.abs(mu)>1-1e-8: y = mu**.5*x
        else: y = Axbar + x; y /= np.sqrt(np.dot(y.conj(),y))

        # constructing a unitary V with first column = y
        V = np.random.random((n,n))+1j*np.random.random((n,n))
        V[:,0] = y
        V,__ = np.linalg.qr(V)

        # update A
        A = np.tensordot(np.tensordot(V.conj().T,A,(1,0)),V.conj(),(1,0))[1:,1:]

        # update U
        Vi = np.eye(N,dtype=complex)
        Vi[N-n:,N-n:] = V
        U = np.tensordot(U,Vi,(1,0))

    return U

a=np.array([0,1,1,0])
a.shape = (2,2)
print a

u= takagi_for_unitary(a)
print u
