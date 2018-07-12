from sympy import *
from sympy.abc import alpha, beta, omega, phi, theta, xi, zeta
from sympy.physics.secondquant import Commutator as Com

class Commutators():
    
    def __init__(self, bmodes, amodes, modetransform):

        self.amodes=amodes
        self.bmodes=bmodes
        self.n=len(bmodes[:,0])
        self.modetrans=modetransform
        
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
    def c(self,A,B):
        print('\nComutator is')
        print('[', end='')
        pprint(A)
        print(',', flush=True)
        pprint(B)
        print(']', flush=True)
        print('\n Which evaluates to:')
        pprint(Com(A,B))

        return Com(A,B)

    def constructmodeops(self):
        return self.d
