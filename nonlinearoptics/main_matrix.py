from sympy import Matrix, init_printing, pi, pprint, symbols, simplify, expand, factor, apart, cancel
from sympy.abc import alpha, beta, omega, phi, zeta
from sympy.physics.secondquant import B, Dagger, BKet, NO
from sympy import *
from makeelements import Makeelements
from makeelements import *
from operatoralg import Operatoralg
import mpmath as mp
import numpy
import math
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
################################################################
##### open data file

g4file = open("g4data.txt", "w+")

g4file.write('xi phi g4\n')

# set increments
xi_start = 0.0
xi_end = 1.0
xi_incr = 0.5

xi_steps = math.ceil((xi_end - xi_start) / float(xi_incr))

phi_start = 0.0
phi_end = float(2 * pi)
phi_incr = float(0.1 * pi)

phi_steps = math.ceil((phi_end - phi_start) / float(phi_incr))

xi_val = xi_start
phi_val = phi_start
####################################

for i in range(0, xi_steps):

    # reset phi for each xi
    phi_val = phi_start

    f_g4norm = g4norm

    #pprint(f_g4norm)
    #f_g4norm = factor(g4norm)

    #print(f_g4norm.free_symbols)
    thing = [None] * len(m.xi)

    # to subs xi
    #xivalue = 0.0

    for i in range(0, len(m.xi)):
        thing[i] = f_g4norm.subs(m.xi[i], xi_val)
        #pprint(thing[i])
    final_g4norm = thing[(len(m.xi) - 1)]

    #pprint(final_g4norm)

    #to subs phi

    # pass phi as arg now
    func = lambdify(m.phi, final_g4norm, "numpy")

    for j in range(0, phi_steps):
        #print('xi= %f, phi= %f, g4= %f' % (xi_val, phi_val, func(phi_val)))
        g4file.write('%f %f %f\n' % (xi_val, phi_val, func(phi_val)))
        phi_val = phi_val + phi_incr

    xi_val = xi_val + xi_incr

g4file.close()
# close data file

print('ran for %d xi vals and %d phi vals for each xi, %d total points' %
      (xi_steps, phi_steps, xi_steps * phi_steps))
exit()


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
#f = 14

H = Matrix([[0, 0, 0, f], [0, 0, (f), 0], [0, conjugate(f), 0, 0],
            [conjugate(f), 0, 0, 0]])

pprint(H)

print('Grouping signal and idler terms,')

A = MatrixSymbol('A', 2, 2)

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
        #Ablk = symbols('A')
        arg = (-I * BlockMatrix([[eye(2), zeros(2)], [zeros(2), eye(2)]]))
        arg = (-I * diag(eye(1), -eye(1)) * Matrix([[0, A], [A, 0]]))
        Hmat = Matrix(H)
        arg = (-I * diag(eye(2), -eye(2)) * Hmat)
    print('arg=')
    pprint(arg)
    return exp(arg)


#print('exp(-iKH')
#exph = (expik(H))
#pprint(exph)
#pprint(simplify(exph))

#
#
#
#

pprint(Ha)

print('type H %s, type Ha %s' % (type(H), type(Ha)))

expha = expik(Ha)
pprint(expha)
simpexpha = simplify(factor(expha))
pprint(simpexpha)


#
#
def cpp_generator(mat, name, label, args=1):
    update_cpp = '\nvoid {0}::{1}_update()'.format(label, name) + '{'
    for n in range(mat.shape[1]):
        for m in range(mat.shape[0]):
            expr = (mat[m, n])
            #symbs = expr.free_symbols
            c = ccode(expr)  #, dereference=args)
            update_cpp += '\n_{0}({1}, {2}) = {3};'.format(name, m, n, c)
    update_cpp += '\n};'
    return update_cpp


simp = simpexpha.subs(f, 13)
#simp = simpexpha
pprint(simp)

label = 'M_matrix'
mat = simp
name = 'F'
code = cpp_generator(mat, name, label)
print(code)


def matrixprinter(outfile, mat):
    nmat = mp.matrix(N(mat))
    outfile.write(str(mat.shape[1]) + ' ')
    outfile.write(str(mat.shape[0]) + '\n')
    for j in range(0, mat.shape[1]):
        for i in range(0, mat.shape[0]):

            print(nmat[i, j])
            outfile.write(str(nmat[i, j]) + ' ')
        outfile.write('\n')
    return 0


with open('temp.txt', 'w') as f:
    #f.write(cxxcode(N(simp)))
    #matrixprinter(f, simp)
    mpmatout = (mp.matrix(N(simp)))
    pprint(mpmatout)
    matrixprinter(f, simp)
"""
f_g4norm = factor(g4norm)

thing = [None] * len(m.xi)

# to subs xi
xivalue = 0.0

for i in range(0, len(m.xi)):
    thing[i] = f_g4norm.subs(m.xi[i], xivalue)

final_g4norm = thing[(len(m.xi) - 1)]

pprint(final_g4norm)

#to subs phi

# pass phi as arg now
func = lambdify(m.phi, final_g4norm)

print(func(0))
"""
#

#
#
#
#
#
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
