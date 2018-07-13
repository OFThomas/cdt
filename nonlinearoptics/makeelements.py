from sympy import (BlockMatrix, I, Matrix, MatrixSymbol, conjugate, cos, cosh,
                   exp, eye, pprint, relational, sin, sinh, sqrt, symbols,
                   zeros)
from sympy.abc import xi
from sympy.physics.quantum import TensorProduct
from sympy.physics.secondquant import Dagger


class Makeelements():
    def __init__(self, nspace, nspectral, modes_in, sq1_modes, sq2_modes,
                 phasemode, phaseangle, bsmodes, bsangle):

        self.nspace = nspace
        self.n = nspace * nspectral
        self.nspectral = nspectral
        self.a = modes_in

        # phase shift
        omega = [None] * (self.nspectral) * 2
        for i in range(0, nspectral):
            omega[i] = symbols('omega%d' % (i))

    # omega[0]=pi/2

    # ########## Define squeezing symbols
        xi = [None] * (self.n**2)
        for i in range(0, self.n**2):
            xi[i] = symbols('xi%d' % (i))

        # modes 0 & 1
        xi[1] = 1
        xi[2] = 2

        # modes 2 & 3
        xi[5] = 1
        xi[6] = 2

        # make active s block anti-diag
        xi[0] = 0
        xi[3] = 0
        xi[4] = 0
        xi[7] = 0

        self.modes = self.makemodes(self.a)

        self.block = self.makeblock()

    def makeblock(self):
        # block form for symplectic matrix
        self.al = MatrixSymbol('alpha', self.n, self.n)
        self.be = MatrixSymbol('beta', self.n, self.n)

        # make block form for symplectic matrix
        self.block = BlockMatrix([[self.al, self.be],
                                  [conjugate(self.be),
                                   conjugate(self.al)]])
        return self.block

    def makemodes(self, symb):
        # Make optical modes
        ablk = MatrixSymbol('a', self.n, 1)
        ablkdag = MatrixSymbol('adag', self.n, 1)
        # build block matrix for modes
        bmodes = BlockMatrix([[ablk], [ablkdag]])
        # make a size n vector for a
        modematrix = Matrix(self.n, 1, lambda i, j: symb[i])
        dagmodematrix = Matrix(self.n, 1, lambda i, j: Dagger(symb[i]))
        # substitution for a modes
        modes = self.matsubs(bmodes, ablk, modematrix, ablkdag, dagmodematrix,
                             0, 0)
        return modes

    def makeps(self, mode1, phaseangle):
        # Phase shifter
        m1 = mode1 * self.nspectral
        phasespace = Matrix(
            self.n, self.n,
            lambda i, j: exp(I * phaseangle[i % self.nspectral]) if ((i == j) and ((m1 <= i) and (i < m1 + self.nspectral))) else 1 if (i == j) else 0
        )
        # symplectic phase shifer
        mps = Matrix(
            self.matsubs(self.block, self.al, phasespace, self.be,
                         zeros(self.n), 0, 0))
        return mps

    def makebs(self, mode, theta):
        # Make beamsplitter
        beamspace = Matrix(
            self.nspace, self.nspace,
            lambda i, j: cos(theta) if ((i == j) and ((i == mode[0]) or (i == mode[1]))) else sqrt(-1) * sin(theta) if (i + j) == (mode[0] + mode[1]) and ((i == mode[0]) or (i == mode[1])) else 1 if (i == j) and ((i != mode[0]) or (i != mode[1])) else 0
        )
        beamsplitter = TensorProduct(beamspace, eye(self.nspectral))
        # symplectic beamsplitter
        mbs = Matrix(
            self.matsubs(self.block, self.al, beamsplitter, self.be,
                         zeros(self.n), 0, 0))
        return mbs

    def makesinglesq(self):
        # Single mode squeezer
        c1 = Matrix(self.n, self.n,
                    lambda i, j: c(xi[i % self.nspectral]) if i == j else 0)
        s1 = Matrix(self.n, self.n,
                    lambda i, j: s(xi[i % self.nspectral]) if i == j else 0)
        return c1, s1

    def makesq(self, mode, sqparam):
        # Two mode squeezer
        # for off-diagonal
        n = self.n
        nspec = self.nspectral
        # for diagonal
        m1s = mode[0] * nspec
        m2s = mode[1] * nspec
        # diagonal part cosh on modes, I elsewhere
        c2 = Matrix(
            n, n,
            lambda i, j: cosh(sqparam[i + m2s]) if (i == j) and ((m1s <= i < m1s + nspec) or (m2s <= i < m2s + nspec)) else 1 if i == j else 0
        )
        # off diag sinh on modes, zero else
        s2 = Matrix(
            n, n,
            lambda i, j: sinh(sqparam[2 * i + j - m1s]) if (((m1s <= i < m1s + nspec) or (m2s <= i < m2s + nspec)) and ((m1s <= j < m1s + nspec) or (m2s <= j < m2s + nspec))) and (abs(i - j) == nspec) and (i != j) else 0
        )
        # symplectic sqs
        ms01 = Matrix(self.matsubs(self.block, self.al, c2, self.be, s2, 0, 0))
        return ms01

    def matsubs(self, mat, aold, anew, bold, bnew, cold, cnew):
        temp1 = mat.subs(aold, anew)
        temp2 = temp1.subs(bold, bnew)
        temp3 = temp2.subs(cold, cnew)
        return temp3

    def justdoitplease(self, transform, modes, showmodes):
        modetransform = transform * modes
        pprint(
            relational.Eq(
                Matrix(modes[0:showmodes]),
                Matrix(modetransform[0:showmodes])))
        return modetransform
