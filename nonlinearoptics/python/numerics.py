from sympy import *
import mpmath as mp
import numpy as np

class Matrixops():
    def __init__(self, matrix):
        self.matrix = matrix
        self.diag = self.svd()

    def svd(self):

        print('\n Numerics \n')

        pprint(m.sq1)
        S1 = m.sq1[range(sq1_mode1, sq1_mode2 + 1),
                   range(sq1_mode1 + n, sq1_mode2 + n + 1)]

        pprint(S1)
        smat = Matrix(N(S1))

        smat = mp.matrix(2 * n)
        imax, jmax = transform.shape

        for j in range(0, jmax):
            for i in range(0, imax):
                smat[i, j] = mp.mpmathify(N(transform[i, j]))

        # smat=Matrix(N(msq1))
        print('\n Squeeze matrix\n')
        # pprint(smat)

        # ######################### svd ##########################
        print('\nnumerical Svd')
        ufloat, sfloat, vfloat = mp.svd_c(smat)

        stemp = mp.matrix(2 * n)
        for i in range(0, 2 * n):
            stemp[i, i] = sfloat[i]

        vtemp = mp.matrix(2 * n)
        for j in range(0, 2 * n):
            for i in range(0, 2 * n):
                vtemp[i, j] = vfloat[i, j]

        # set numbers <10^-15 = 0
        u = mp.chop(ufloat)
        s = mp.chop(stemp)
        v = mp.chop(vtemp)

        # reconstruct to compare
        sqreconstruct = mp.chop(u * s * v)

        print('\nAfter svd doing u*s*v\n')
        # pprint(sqreconstruct)

        print('\nOriginal sq matrix\n')
        # pprint(smat)

        diff = sqreconstruct - smat
        # pprint(mp.chop(diff))

        print(mp.norm(diff, p='inf'))

        if mp.norm(diff, p='inf') <= 10**-5:
            print('\n#############################################')
            print('#######        Close enough     #############')
            print('#############################################\n')

    def makematrixexact(self, matrix):
        imax, jmax = matrix.shape
        newmat = matrix[:, :]
        for j in range(0, jmax):
            for i in range(0, imax):
                newmat[i, j] = nsimplify(matrix[i, j], full=True)

        return matrix

    def takagi_for_unitary(self, A):
        ### takes a unitary matrix A such that A^T = A, ###
        ###    returns unitary U such that A = U U^T    ###
        N = A.shape[0]
        U = np.eye(N)

        while A.shape[0] > 0:
            n = A.shape[0]
            # constructing a vector v such that A*ybar = y
            x = np.random.rand(n) + 1j * np.random.rand(n)
            x /= np.sqrt(np.dot(x.conj(), x))
            Axbar = np.tensordot(A, x.conj(), (1, 0))
            mu = np.dot(x.conj(), Axbar)

            if np.abs(mu) > 1 - 1e-8:
                y = mu**.5 * x
            else:
                y = Axbar + x
                y /= np.sqrt(np.dot(y.conj(), y))

            # constructing a unitary V with first column = y
            V = np.random.random((n, n)) + 1j * np.random.random((n, n))
            V[:, 0] = y
            V, __ = np.linalg.qr(V)
            # update A
            A = np.tensordot(
                np.tensordot(V.conj().T, A, (1, 0)), V.conj(), (1, 0))[1:, 1:]
            # update U
            Vi = np.eye(N, dtype=complex)
            Vi[N - n:, N - n:] = V
            U = np.tensordot(U, Vi, (1, 0))

        return U
