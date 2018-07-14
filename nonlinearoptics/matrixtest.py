from sympy import Matrix, init_printing, pi, pprint, symbols, simplify
from sympy.abc import alpha, beta, omega, phi, zeta
from sympy.physics.secondquant import B, Dagger, BKet, NO

from makeelements import Makeelements
from operatoralg import Commutators

init_printing(use_unicode=True)


def printc(matrix):
    print()
    pprint(matrix[0:n, 0:n])


def prints(matrix):
    print()
    pprint(matrix[0:n, n:2 * n])


# ################################# start of program

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

# input and output modes
a = [None] * n
b = [None] * n
for i in range(0, nspace):
    for j in range(0, nspectral):
        # B is boson annihilation op
        a[i * (nspectral) + j] = B('%d+%d' % (i * nspectral, j))
        b[i * (nspectral) + j] = B('%d+%d' % (i * nspectral, j))

# ################# Make stuff happen!#######################

# makes a bs, ps, sq1, sq2
m = Makeelements(nspace, nspectral, a, sq1_modes, sq2_modes, phasemode,
                 phaseangle, bsmodes, bsangle)

# ############## numerics ###################################

# phase shift
theta = [None] * (nspectral * 2)
for i in range(0, nspectral):
    theta[i] = symbols('theta%d' % (i), real=True)

xi = [None] * (n**2)
for i in range(0, n**2):
    xi[i] = symbols('r%d' % (i), real=True)
    xi[i] = symbols('r', real=True)

phi = symbols('phi', real=True)

transform = [None] * 4

squeezer1 = m.makesq(mode=[0, 1], sqparam=xi)
squeezer2 = m.makesq(mode=[2, 3], sqparam=xi)
phaseshift = m.makeps(mode1=1, phaseangle=theta)
beamsplitter = m.makebs(mode=[1, 2], theta=phi)

transform[0] = squeezer1
transform[1] = squeezer2 * transform[0]
transform[2] = phaseshift * transform[1]
transform[3] = beamsplitter * transform[2]
"""
pprint(squeezer1)
pprint(squeezer2)
pprint(phaseshift)
pprint(beamsplitter)
"""

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

bmodes = Matrix(m.makemodes(b))
amodes = m.modes[:, 0]

# for characterising
fulltransform = Matrix(m.makeblock())
fulltransform = transform[len(transform) - 1]

modetrans = m.justdoitplease(fulltransform, m.modes, showmodes=n)
# end of characterising

# from operatoralg import commutation stuff
com = Commutators(bmodes, amodes, modetrans, fulltransform)

# use relational to make d final mode transforms
d = com.constructmodeops()

print('\n\n\n')


def correlationfn(indices, modes):
    g = 1
    for i in range(0, len(indices)):
        print('i,j', i, j)
        g = g * modes[i]
    for i in range(len(indices) - 1, -1, -1):
        g = g * modes[n + i]
    return g


c_d0d1 = com.c(d[0], d[1])
c_d0dagd0 = com.c(d[0], d[n + 0])

# mode transforms
pprint(d)

# commutation relations using matrix elements

#pprint(com.matrixel([0, 1], d))
#pprint(com.matrixel([1, 0], d))

#pprint(com.matrixel([n + 0, n + 1], d))
#pprint(com.matrixel([n + 1, n + 0], d))

#pprint(com.matrixel([0, n + 0], d))
pprint(com.matrixel([n + 0, n + 1], d) * com.matrixel([n + 2, n + 3], d))

g4 = 0
counter = 0
print('G4 calc')

for j in range(0, 2 * n):
    for i in range(j + 1, 2 * n):
        print('#############################################')
        print('step i,j', i, j)
        g4 = g4 + com.matrixel([j, i], d)
        pprint(g4)
        counter += 1

print('Final g4 =')
pprint(g4)
print('counter =', counter)

print('simplify')
pprint(simplify(g4))

# ### for real

# def picktwonumbers(listin):

g4 = 0
multi = 1
counter = 0
print('G4 calc')

# make list 0 to 2*n-1
nums = [i for i in range(0, 2 * n)]

for j in range(0, n):
    print('#############################################')
    print('step j', j)

    for i in range(0, len(nums)):
        #       multi=multi*
        1

    g4 = g4 + com.matrixel([j, i], d)
    pprint(g4)
    counter += 1

print('Final g4 =')
pprint(g4)
print('counter =', counter)

print('simplify')
pprint(simplify(g4))
#
