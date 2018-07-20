from sympy import Matrix, init_printing, pi, pprint, symbols, simplify, expand, factor, apart, cancel
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
    xi[i] = symbols('r%d' % (i))
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
if showexplicit == 1:
    fulltransform = Matrix(m.makeblock())
else:
    fulltransform = transform[len(transform) - 1]

modetrans, rel = m.justdoitplease(fulltransform, m.modes, showmodes=n)
# end of characterising

print('rel')
pprint(rel)
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
#pprint(com.matrixel([n + 0, n + 1], d) * com.matrixel([n + 2, n + 3], d))


def amp(a):
    return abs(a)**2


g4 = 0

gam10 = com.ABT(1, 0)
gam32 = com.ABT(3, 2)
gam21 = com.ABT(2, 1)
gam30 = com.ABT(3, 0)
gam20 = com.ABT(2, 0)
gam31 = com.ABT(3, 1)

gam = amp(gam10 * gam32 + gam21 * gam30 + gam20 * gam31)

print('\n gamma \n')
pprint(gam)

bbdag00 = com.BBD(0, 0)
bbdag11 = com.BBD(1, 1)
bbdag22 = com.BBD(2, 2)
bbdag33 = com.BBD(3, 3)

bdiag = bbdag00 * bbdag11 * bbdag22 * bbdag33
print('\nBBdag diag terms\n')
pprint(bdiag)

term=[None]*6

term[0] = amp(gam10) * bbdag33 * bbdag22
term[1] = amp(gam21) * bbdag00 * bbdag33
term[2] = amp(gam20) * bbdag11 * bbdag33
term[3] = amp(gam32) * bbdag00 * bbdag11
term[4] = amp(gam31) * bbdag00 * bbdag22
term[5] = amp(gam30) * bbdag11 * bbdag22

g4 = gam + bdiag + term[0] + term[1] + term[2] + term[3] + term[4] + term[5]

G4 = g4  #.subs(phi,pi/4)
print('\nG4\n')
pprint(G4)
print('\nSimplified\n')
pprint((G4))

print('\n no beamsplitter\n')
G4_no_bs = (g4.subs(phi, 0))
pprint(G4_no_bs)

print('answer')
print('\n\n\n\n help?')

g4norm=(simplify(G4 / G4_no_bs))

pprint((g4norm))

print('\n g4')


"""
d_nobs=d.subs(phi,0)

pprint(d_nobs)
for i in range(0,2*n):
    print()
    pprint(d_nobs[i])
"""


"""
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
