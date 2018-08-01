from sympy import Matrix, init_printing, pi, pprint, symbols, simplify, expand, factor, apart, cancel
from sympy.abc import alpha, beta, omega, phi, zeta
from sympy.physics.secondquant import B, Dagger, BKet, NO
from sympy import *
from makeelements import Makeelements
from makeelements import *
from operatoralg import Operatoralg

init_printing(use_unicode=True)


def printc(matrix):
    print()
    mat = (matrix[0:n, 0:n])
    pprint(mat)
    return mat


def prints(matrix):
    print()
    mat = (matrix[0:n, n:2 * n])
    pprint(mat)
    return mat


# ################################# start of program
showexplicit = 0

# spatial dim
nspace = 4
# spectral dim
nspectral = 1

# make total dimension
n = nspace * nspectral

# ################ Specify modes for unitaries
# 2 mode squeezer on modes 0 & 1
sq1_modes = [0, 1]

# 2 mode sq on modes 2 & 3
sq2_modes = [2, 3]

# Phase shifter
phasemode = 2
phaseangle = pi / 2

# beamsplitter spatial modes
bsmodes = [1, 2]
bsangle = pi / 4

# ########## Define optical modes
print('Spatial modes =', nspace)
print('Spectral modes =', nspectral)
print('Total dim =', n)
print('\n Two mode squeezing on modes, ', sq1_modes[0], ',', sq1_modes[1])
print('Two mode squeezer on modes, ', sq2_modes[0], ',', sq2_modes[1])
print('Phase shift on modes. ', phasemode, 'phase angle', phaseangle)
print('Beamsplitter on modes, ', bsmodes[0], ',', bsmodes[1], 'angle= ',
      bsangle)

# ################# Make stuff happen!#######################

# makes a bs, ps, sq1, sq2
m = Makeelements(nspace, nspectral, sq1_modes, sq2_modes, phasemode,
                 phaseangle, bsmodes, bsangle)

# ############## numerics ###################################

transform = [None] * 4

squeezer1 = m.makesq(mode=[0, 1])  #, sqparam=xi)
squeezer2 = m.makesq(mode=[2, 3])  #, sqparam=xi)
phaseshift = m.makeps(mode1=1)  #, phaseangle=theta)
beamsplitter = m.makebs(mode=[1, 2])  #, theta=phi)

transform[0] = squeezer1
transform[1] = squeezer2 * transform[0]
transform[2] = phaseshift * transform[1]
transform[3] = beamsplitter * transform[2]

# do mode transformation
for i in range(0, len(transform)):
    print('\n mode transform')
    modetrans = m.justdoitplease(transform[i], m.modes, showmodes=n)
    print('\nC matrix')
    printc(transform[i])
    print('\nS matrix')
    prints(transform[i])
    print()

# ################## print s and c matrix
print('\nC matrix')
printc(transform[len(transform) - 1])
print('\nS matrix')
prints(transform[len(transform) - 1])
print()

bmodes = Matrix(m.makemodes(m.b))
amodes = m.modes[:, 0]

# for characterising
if showexplicit == 1:
    fulltransform = Matrix(m.makeblock())
else:
    fulltransform = transform[len(transform) - 1]

modetrans, rel = m.justdoitplease(fulltransform, m.modes, showmodes=n)
# end of characterising

print('rel')
pprint(rel)
# from operatoralg import commutation stuff
opalg = Operatoralg(bmodes, amodes, modetrans, fulltransform)

# use relational to make d final mode transforms
d = opalg.constructmodeops()

print('\n\n\n')

# c_d0d1 = opalg.c(d[0], d[1])
# c_d0dagd0 = opalg.c(d[0], d[n + 0])

# mode transforms
pprint(d)

# commutation relations using matrix elements

print('\n g4')
g4 = opalg.calcg4()
g4_nobs = g4.subs(m.phi, 0)
g4norm = factor(g4) / factor(g4_nobs)

print('Factorised G4 normalised')
pprint(factor(g4norm))


def jordanblock(i, j):
    if i == j - 1:
        return 1
    elif i == j:
        return x[i]
    else:
        return 0


def matrixcosh(x):
    return (exp(x) + exp(-x)) / 2


def matrixsinh(x):
    return (exp(x) - exp(-x)) / 2


def dag(x):
    return conjugate(Transpose(x))


x = symbols('x:10')

print('H')
f = symbols('F')
H = Matrix([[0, 0, 0, f], [0, 0, (f), 0], [0, conjugate(f), 0, 0],
            [conjugate(f), 0, 0, 0]])

pprint(H)

print('Grouping signal and idler terms,')

#A = MatrixSymbol('A', 2, 2)

A = Matrix([[0, f], [(f), 0]])
Ha = BlockMatrix([[zeros(2), A], [(conjugate(A)), zeros(2)]])

pprint(H)


def expik(H):
    if (str(type(H)) == "<class 'sympy.matrices.dense.MutableDenseMatrix'>"):
        print('matrix')
        arg = (-I * diag(eye(2), -eye(2)) * H)
    elif (str(type(H)) ==
          "<class 'sympy.matrices.expressions.blockmatrix.BlockMatrix'>"):
        print('block matrix')
        A = symbols('A')
        arg = (-I * BlockMatrix([[eye(2), zeros(2)], [zeros(2), eye(2)]]))
        arg = (-I * diag(eye(1), -eye(1)) * Matrix([[0, A], [A, 0]]))
    print('arg=')
    pprint(arg)
    return exp(arg)


print('exp(-iKH')
exph = (expik(H))
pprint(exph)
pprint(simplify(exph))

#
#
#
#

pprint(Ha)

print('type H %s, type Ha %s' % (type(H), type(Ha)))

expha = expik(Ha)
pprint(expha)
pprint(simplify(factor(expha)))
"""
matx = Matrix([[x[0], x[1]], [0, x[3]]])

matx = Matrix(3, 3, jordanblock)
pprint(matx)
# ###############################
coshx = matrixcosh(matx)
print('\n Cosh of matrix')
print('\n unsimplified')
# pprint(coshx)

print('\nsimplified')
pprint(simplify(coshx))

# ###############################################
print('\nSinh of matrix')
sinhx = matrixsinh(matx)
print('\nunsimplified')
# pprint(sinhx)

print('\nsimplified')
pprint(simplify(sinhx))

beta=prints(fulltransform)
betadag=conjugate(prints(fulltransform))
print('BB^ =')
pprint(beta*betadag.T)


"""
#
#
#

#
