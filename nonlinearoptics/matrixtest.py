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
"""
# phase shift
theta = [None] * (nspectral * 2)
for i in range(0, nspectral):
    theta[i] = symbols('theta%d' % (i), real=True)

# squeeze param
xi = [None] * (n**2)
for i in range(0, n**2):
    xi[i] = symbols('r%d' % (i))
    xi[i] = symbols('r', real=True)

# beamsplitter angle
phi = symbols('phi', real=True)
"""
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
H = Matrix([[0, 0, 0, f], [0, 0, Transpose(f), 0], [0, conjugate(f), 0, 0],
            [dag(f), 0, 0, 0]])

pprint(H)

print('Grouping signal and idler terms,')
#H=Matrix([

#A = MatrixSymbol('A', 2, 2)

A = Matrix([[0, f], [Transpose(f), 0]])
Ha = BlockMatrix([[0, A], [(conjugate(A)), 0]])

pprint(H)


def expik(H):
    arg = (-I * diag(1, 1, -1, -1) * H)
    print('arg=')
    pprint(arg)
    return exp(arg)


print('exp(-iKH')
exph = (expik(H))
pprint(exph)
pprint(simplify(exph))
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

#pprint(g4.subs(m.phi, 0))

d_nobs=d.subs(phi,0)

pprint(d_nobs)
for i in range(0,2*n):
    print()
    pprint(d_nobs[i])


cresult=[None]*12

print('\nd0dag')
cresult[0]=(com.matrixel([0 + n, 1 + n], d))
cresult[1]=(com.matrixel([0 + n, 2 + n], d))
cresult[2]=(com.matrixel([0 + n, 3 + n], d))

cresult[3]=(com.matrixel([0 + n, 3], d))
cresult[4]=(com.matrixel([0 + n, 2], d))
cresult[5]=(com.matrixel([0 + n, 1], d))

print('\nd1dag')
cresult[6]=(com.matrixel([1 + n, 2 + n], d))
cresult[7]=(com.matrixel([1 + n, 3 + n], d))

cresult[8]=(com.matrixel([1 + n, 3], d))
cresult[9]=(com.matrixel([1 + n, 2], d))

print('\nd2dag')
cresult[10]=(com.matrixel([2 + n, 3 + n], d))
cresult[11]=(com.matrixel([2 + n, 3], d))

print(cresult)

g4expr=1
for i in range(n-1,-1,-1):
    print('\nd %d dagd %d\n' % (i,i))
    ddagd=expand(d[i+n]*d[i])
    pprint(expand(ddagd))
    g4expr=g4expr*ddagd

print('\n\n g4 expr=')
pprint((g4expr))
"""

#
