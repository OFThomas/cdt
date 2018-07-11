from sympy import *
from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum import TensorProduct 
from mpmath import mp
from sympy.abc import alpha, beta, xi, zeta

from makeelements import *

init_printing(use_unicode=True)

################################## start of program

#spatial dim
nspace=4

#spectral dim 
nspectral=1

#make total dimension
n=nspace*nspectral

################# Specify modes for unitaries
# 2 mode squeezer on modes 0 & 1
sq1_mode1,sq1_mode2= 0,1

# 2 mode sq on modes 2 & 3
sq2_mode1,sq2_mode2=2,3

# Phase shifter
phasemode=2
phaseangle=pi/2

# beamsplitter spatial modes
bsmode1, bsmode2= 1,2
bsangle=pi/4

#cosh & sinh placeholders
c,s = symbols('c s')

########### Define optical modes
print('Spatial modes =', nspace)
print('Spectral modes =', nspectral)
print('Total dim =', n)
print('\n Two mode squeezing on modes, ', sq1_mode1,',', sq1_mode2)
print('Two mode squeezer on modes, ', sq2_mode1, ',', sq2_mode2)
print('Phase shift on modes. ', phasemode, 'phase angle', phaseangle)
print('Beamsplitter on modes, ', bsmode1, ',', bsmode2, 'angle= ', bsangle)


################## Make stuff happen!#######################

# makes a bs, ps, sq1, sq2
m = Makeelements(n,nspectral, sq1_mode1, sq1_mode2,
                    sq2_mode1, sq2_mode2, phasemode, phaseangle,
                    bsmode1, bsmode2, bsangle)
 

#do mode transformation
transform=m.bs*m.ps*m.sq2*m.sq1 
m.justdoitplease(transform,m.modes)

#############################################################
############### numerics ###################################
print('\n Numerics \n')

pprint(m.sq1)
S1=m.sq1[range(sq1_mode1,sq1_mode2+1),range(sq1_mode1+n,sq1_mode2+n+1)]

pprint(S1)
smat=Matrix(N(S1))

smat=mp.matrix(2*n)
imax,jmax=transform.shape

for j in range(0,jmax):
    for i in range(0,imax):
        smat[i,j]=mp.mpmathify(N(transform[i,j]))


#smat=Matrix(N(msq1))
print('\n Squeeze matrix\n')
#pprint(smat)

########################## svd ##########################
print('\nnumerical Svd')
ufloat,sfloat,vfloat = mp.svd_c(smat)

stemp=mp.matrix(2*n)
for i in range(0,2*n):
    stemp[i,i]=sfloat[i]

vtemp=mp.matrix(2*n)
for j in range(0,2*n):
    for i in range(0,2*n):
        vtemp[i,j]=vfloat[i,j]

# set numbers <10^-15 = 0 
u=mp.chop(ufloat)
s=mp.chop(stemp)
v=mp.chop(vtemp)

# reconstruct to compare 
sqreconstruct=mp.chop(u*s*v)

print('\nAfter svd doing u*s*v\n')
#pprint(sqreconstruct)

print('\nOriginal sq matrix\n')
#pprint(smat)

diff=sqreconstruct-smat 
#pprint(mp.chop(diff))

print(mp.norm(diff,p='inf'))

if (mp.norm(diff,p='inf')<=10**-5):
    print('\n#############################################')
    print('#######        Close enough     #############')
    print('#############################################\n')

# 
# pprint(m.makeps(0,[2],2))
# bs
# pprint(m.makebs(0,1,2))


