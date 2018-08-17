from sympy import *
from sympy.abc import alpha, beta, omega, theta, xi, zeta
from sympy.physics.secondquant import NO
from sympy.physics.secondquant import Commutator as Com
from sympy.physics.secondquant import wicks

from makeelements import *


class Operatoralg():
    def __init__(self, nspec,bmodes, amodes, modetransform, fulltransform):
        self.fulltransform = fulltransform
        self.amodes = amodes
        self.bmodes = bmodes
        self.nspec=nspec
        #        self.n = len(bmodes[:, 0])
        self.n = int((amodes.shape)[0])
        print('n', self.n)

        self.modetrans = modetransform

        self.moderel = [None] * self.n
        # make new modes b = Ma
        #for i in range(0, self.n):
        #pprint(type(self.modetrans))
        #self.moderel[i] = Eq(self.bmodes[i, 0], self.modetrans[i, 0])

        # whole mode transform matrix
        #pprint(Eq(self.bmodes[:, 0], self.modetrans[:, 0]))

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

    def ABT(self, i0, i1, transform):
        ft = transform
        n = int(0.5 * self.n)
        if (i0 < n) and (i1 < n):
            a = 0
            b = n
        elif (i0 >= n) and (i1 >= n):
            a = n
            b = 2 * n
            i0 = i0 - n
            i1 = i1 - n

        matel = ft[a:b, 0:n] * ft[a:b, n:2 * n].T
        res=0
        for i in range(0,self.nspec):
            res=res+matel[(i0*self.nspec)+i,(i1*self.nspec)+i]
        return res

    def BBD(self, i0, i1, transform):
        ft = transform
        n = int(0.5 * self.n)
        matel = ft[0:n, n:2 * n] * ft[n:2 * n, 0:n].T
        res=0
        for i in range(0,self.nspec):
            res=res+matel[(i0*self.nspec)+i,(i1*self.nspec)+i]
        return res

    def matrixel(self, index, A):
        n = int(0.5 * self.n)
        fulltransform = self.fulltransform
        print('\ncom d[', index[0], '], d[', index[1], ']')
        self.c(A[index[0]], A[index[1]])
        print()
        # alpha*beta.T
        if (index[0] < n) and (index[1] < n) or ((index[0] >= n) and
                                                 (index[1] >= n)):

            #print('Matrix element')
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

    def calcg4(self, transform):
        def amp(a):
            return abs(a)**2
        
        trans=transform
        g4 = 0

        gam10 = self.ABT(1, 0,trans)
        gam32 = self.ABT(3, 2,trans)
        gam21 = self.ABT(2, 1, trans)
        gam30 = self.ABT(3, 0,trans)
        gam20 = self.ABT(2, 0, trans)
        gam31 = self.ABT(3, 1,trans)

        print('gamma')
        pprint(gam10 * gam32)
        print('gamma')
        pprint(gam21 * gam30)
        print('gamma')
        pprint(gam20 * gam31)
        gam = amp(gam10 * gam32 + gam21 * gam30 + gam20 * gam31)

        bbdag00 = self.BBD(0, 0, trans)
        bbdag11 = self.BBD(1, 1, trans)
        bbdag22 = self.BBD(2, 2, trans)
        bbdag33 = self.BBD(3, 3, trans)
        print('beta beta 00', self.BBD(0,0, trans))
        bdiag = bbdag00 * bbdag11 * bbdag22 * bbdag33

        term = [None] * 6

        term[0] = amp(gam10) * bbdag33 * bbdag22
        term[1] = amp(gam21) * bbdag00 * bbdag33
        term[2] = amp(gam20) * bbdag11 * bbdag33
        term[3] = amp(gam32) * bbdag00 * bbdag11
        term[4] = amp(gam31) * bbdag00 * bbdag22
        term[5] = amp(gam30) * bbdag11 * bbdag22

        g4 = gam + bdiag + term[0] + term[1] + term[2] + term[3] + term[4] + term[5]
        return g4

        """
        #print('\n no beamsplitter\n')
        G4_no_bs = (g4.subs(phi, 0))
        #pprint(G4_no_bs)

        #print('answer')
        #print('\n\n\n\n help?')

        g4norm = (simplify(g4 / G4_no_bs))
        #pprint((g4norm))
        print('\n gamma \n')
        pprint(gam)
        print('\nBBdag diag terms\n')
        pprint(bdiag)
        print('g4 not normalised')
        pprint(g4)
        pprint(G4_no_bs.subs(phi, 1))
        print('g4calc in opalg')

        print('full transform subs')
        pprint(self.fulltransform)
        print('phi=0')
        pprint(self.fulltransform.subs(phi, 0))
        return g4
        """
