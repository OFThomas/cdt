from mpmath import mp
from sympy import *
from sympy.abc import alpha, beta, omega, phi, theta, xi, zeta
from sympy.physics.secondquant import Commutator as Com

from makeelements import *
from operatoralg import *

init_printing(use_unicode=True)


def printc(matrix):
    print()
    pprint(matrix[0:n, 0:n])


def prints(matrix):
    print()
    pprint(matrix[0:n, n:2 * n])


################################## start of program

#spatial dim
nspace = 4

#spectral dim
nspectral = 8

#make total dimension
n = nspace * nspectral

################# Specify modes for unitaries
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

########### Define optical modes
print('Spatial modes =', nspace)
print('Spectral modes =', nspectral)
print('Total dim =', n)
print('\n Two mode squeezing on modes, ', sq1_modes[0], ',', sq1_modes[1])
print('Two mode squeezer on modes, ', sq2_modes[0], ',', sq2_modes[1])
print('Phase shift on modes. ', phasemode, 'phase angle', phaseangle)
print('Beamsplitter on modes, ', bsmodes[0], ',', bsmodes[1], 'angle= ',
      bsangle)

a = [None] * n
b = [None] * n
for i in range(0, nspace):
    for j in range(0, nspectral):
        #self.a[i*(self.nspectral)+j]=symbols('a%d%d' % (i, j))
        a[i * (nspectral) + j] = B('%d+%d' % (i * nspectral, j))
        b[i * (nspectral) + j] = B('%d+%d' % (i * nspectral, j))
# ################# Make stuff happen!#######################

# makes a bs, ps, sq1, sq2
m = Makeelements(nspace, nspectral, a, sq1_modes, sq2_modes, phasemode,
                 phaseangle, bsmodes, bsangle)

# do mode transformation
#transform = m.bs * m.ps * m.sq2 * m.sq1
#m.justdoitplease(transform, m.modes, showmodes=2 * n)

# ############## numerics ###################################

# phase shift
theta = [None] * (nspectral * 2)
for i in range(0, nspectral):
    theta[i] = symbols('theta%d' % (i))

#theta[0]=symbols('theta')

########### Define squeezing symbols
xi = [None] * (n**2)
for i in range(0, n**2):
    xi[i] = symbols('r%d' % (i))

#modes 0 & 1
#xi[1]=1
#xi[2]=2

#modes 2 & 3
#xi[5]=1
#xi[6]=2

#make active s block anti-diag
#xi[0] = 0
#xi[3] = 0
#xi[4] = 0
#xi[7] = 0

transform = [None] * 4

squeezer1 = m.makesq(mode=[0, 1], sqparam=xi)
squeezer2 = m.makesq(mode=[2, 3], sqparam=xi)
phaseshift = m.makeps(mode1=1, phaseangle=theta)
beamsplitter = m.makebs(mode=[1, 2], theta=phi)

transform[0] = squeezer1
transform[1] = squeezer2 * transform[0]
transform[2] = phaseshift * transform[1]
transform[3] = beamsplitter * transform[2]

#pprint(squeezer1)
#pprint(squeezer2)
#pprint(phaseshift)
#pprint(beamsplitter)

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

com = Commutators(bmodes, amodes, modetrans)

d = com.constructmodeops()

print('\n\n\n')

c_d0d1 = com.c(d[0], d[1])
c_d0dagd0 = com.c(d[0], d[n + 0])

#
