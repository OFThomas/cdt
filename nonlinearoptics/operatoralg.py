from sympy import *
from sympy.abc import alpha, beta, omega, phi, theta, xi, zeta
from sympy.physics.secondquant import NO
from sympy.physics.secondquant import Commutator as Com
from sympy.physics.secondquant import wicks


class Commutators():
    def __init__(self, bmodes, amodes, modetransform, fulltransform):
        self.fulltransform = fulltransform
        self.amodes = amodes
        self.bmodes = bmodes
        self.n = len(bmodes[:, 0])
        self.modetrans = modetransform

        self.moderel = [None] * self.n
        # make new modes b = Ma
        for i in range(0, self.n):
            self.moderel[i] = Eq(self.bmodes[i, 0], self.modetrans[i, 0])

        # whole mode transform matrix
        pprint(Eq(self.bmodes[:, 0], self.modetrans[:, 0]))

        # for formatting
        print('\n\n\n\n\n\n\n\n\n\n\n\n\n')

        # print new b modes
        for i in range(0, self.n):
            pprint(self.moderel[i])

        # for commutators
        self.d = self.modetrans[:, 0]
        #return self.d

    # commutator
    def c(self, A, B):
        print('\nComutator is')
        print('[', end='')
        pprint(A)
        print(',', flush=True)
        pprint(B)
        print(']', flush=True)
        print('\n Which evaluates to:')
        pprint(Com(A, B))

        return Com(A, B)

    def matrixel(self, index, A):
        n = int(0.5 * self.n)
        fulltransform = self.fulltransform
        print('com d[', index[0], '], d[', index[1], ']')
        self.c(A[index[0]], A[index[1]])
        print()
        # alpha*beta.T
        if (index[0] < n) and (index[1] < n) or ((index[0] >= n) and
                                                 (index[1] >= n)):
            print('Matrix element')
            a = 0
            b = n
            # conjugate alpha*beta.T
            if (index[0] >= n) and (index[1] >= n):
                a = n
                b = 2 * n
                index[0] = index[0] - n
                index[1] = index[1] - n

            matelement = fulltransform[a:b, 0:n] * fulltransform[a:b, n:2 *
                                                                 n].T
            return (matelement[index[0], index[1]] -
                    matelement[index[1], index[0]])

        # dagger non-dagger
        if (index[0] >= n) and (index[1] < n) or (index[0] < n) and (index[1]
                                                                     >= n):
            print('Matrix element')
            a = 0
            b = n
            matelement = fulltransform[a:b, 0:n] * fulltransform[n:2 * n, n:2 *
                                                                 n].T
            mat2 = fulltransform[0:n, n:2 * n] * fulltransform[n:2 * n, 0:n].T
            # to fix minus sign for n+0 -n =0
            if index[0] >= n:
                index[0] = index[0] - n
                ans = (
                    matelement[index[0], index[1]] - mat2[index[0], index[1]])
                return (-ans)
            elif index[1] >= n:
                index[1] = index[1] - n
                ans = (
                    matelement[index[0], index[1]] - mat2[index[0], index[1]])
                return (ans)

    def constructmodeops(self):
        return self.d

    def dowicks(self, A):
        pprint(simplify(A))
        return (A)

    def g4(self, mat):
        alpha = mat[0:self.n, 0:self.n]
        beta = mat[0:self.n, self.n:2 * self.n]
        g4 = 0
        #g4=abs()
        return g4
