from sympy import (BlockMatrix, block_collapse, I, Matrix, MatrixSymbol,
                   conjugate, cos, cosh, exp, eye, pprint, relational, sin,
                   sinh, sqrt, symbols, zeros, Eq)
from sympy.abc import xi
from sympy.physics.quantum import TensorProduct
from sympy.physics.secondquant import Dagger, B


class Makeelements():
    def __init__(self, nspace, nspectral, sq1_modes, sq2_modes, phasemode,
                 phaseangle, bsmodes, bsangle):

        self.nspace = nspace
        self.n = nspace * nspectral
        self.nspectral = nspectral

        # phase shift
        omega = [None] * (self.nspectral) * 2
        for i in range(0, nspectral):
            omega[i] = symbols('omega%d' % (i))

        # input and output modes
        self.a = [None] * self.n
        self.b = [None] * self.n
        for i in range(0, self.nspace):
            for j in range(0, self.nspectral):
                # B is boson annihilation op
                self.a[i * (self.nspectral) + j] = B(
                    '%d+%d' % (i * self.nspectral, j))
                self.b[i * (self.nspectral) + j] = B(
                    '%d+%d' % (i * self.nspectral, j))

        # phase shift
        self.theta = [None] * (self.nspectral * 2)
        for i in range(0, self.nspectral):
            self.theta[i] = symbols('theta%d' % (i), real=True)

        # squeeze param
        self.xi = [None] * (self.n**2)
        for i in range(0, self.n**2):
            self.xi[i] = symbols('r%d' % (i))
            self.xi[i] = symbols('r', real=True)

        # beamsplitter angle
        self.phi = symbols('phi', real=True)

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

    def makeps(self, mode1):  #, phaseangle):
        # Phase shifter
        m1 = mode1 * self.nspectral
        phasespace = Matrix(
            self.n, self.n,
            lambda i, j: exp(I * self.theta[i % self.nspectral]) if ((i == j) and ((m1 <= i) and (i < m1 + self.nspectral))) else 1 if (i == j) else 0
        )
        # symplectic phase shifer
        mps = Matrix(
            self.matsubs(self.block, self.al, phasespace, self.be,
                         zeros(self.n), 0, 0))
        return mps

    def makebs(self, mode):  #, theta):
        # Make beamsplitter

        beamspace = Matrix(
            self.nspace, self.nspace,
            lambda i, j: cos(self.phi) if ((i == j) and ((i == mode[0]) or (i == mode[1]))) else sqrt(-1) * sin(self.phi) if (i + j) == (mode[0] + mode[1]) and ((i == mode[0]) or (i == mode[1])) else 1 if (i == j) and ((i != mode[0]) or (i != mode[1])) else 0
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

    def makesq(self, mode):  #, sqparam):
        # Two mode squeezer
        # for off-diagonal
        def diag(i, j):
            if (i == j) and ((m1s <= i < m1s + nspec) or
                             (m2s <= i < m2s + nspec)):
                return cosh(self.xi[i + m2s])
            elif i == j:
                return 1
            else:
                return 0

        def off_diag(i, j):
            if (((m1s <= i < m1s + nspec) or (m2s <= i < m2s + nspec)) and
                ((m1s <= j < m1s + nspec) or (m2s <= j < m2s + nspec))) and (
                    abs(i - j) == nspec) and (i != j):
                return sinh(self.xi[2 * i + j - m1s])
            else:
                return 0

        n = self.n
        nspec = self.nspectral
        # for diagonal
        m1s = mode[0] * nspec
        m2s = mode[1] * nspec

        # diagonal part cosh on modes, I elsewhere
        c2 = Matrix(n, n, diag)
        # off diag sinh on modes, zero else
        s2 = Matrix(n, n, off_diag)
        # symplectic sqs
        ms01 = (self.matsubs(self.block, self.al, c2, self.be, s2, 0, 0))
        return Matrix(ms01)

    def matsubs(self, mat, aold, anew, bold, bnew, cold, cnew):
        temp1 = mat.subs(aold, anew)
        temp2 = temp1.subs(bold, bnew)
        temp3 = temp2.subs(cold, cnew)
        return block_collapse(temp3)

    def justdoitplease(self, transform, modes, showmodes):
        modetransform = transform * modes
        # pprint(relational.Eq(Matrix(modes[0:showmodes]),Matrix(modetransform[0:showmodes])))

        rel = Eq(modes, (transform * modes))
        # pprint(rel)

        # whole mode transform matrix
        # pprint(Eq(modes[:, 0], modetransform[:, 0]))
        return modetransform, rel
