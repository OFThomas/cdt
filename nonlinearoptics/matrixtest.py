from sympy import *
from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum import TensorProduct 
from sympy.abc import alpha, beta, xi, zeta
init_printing(use_unicode=True)

def makebs(mode1,mode2, n, theta):
    #beamsplitter
    passive=Matrix(n,n, lambda i,j: cos(theta) if ((i==j)and((i==mode1)or(i==mode2))) 
            else  sqrt(-1)*sin(theta) if (i+j)==(mode1+mode2)and((i==mode1)or(i==mode2)) 
            else 1 if (i==j)and((i!=mode1)or(i!=mode2)) else 0)
    return passive

def makesinglesq(n, nspec):
    #Single mode squeezer
    c1=Matrix(n,n, lambda i,j: c(xi[i%nspec]) if i==j else 0)
    s1=Matrix(n,n, lambda i,j: s(xi[i%nspec]) if i==j else 0)
    return c1, s1
    
def maketwomodesq(mode1,mode2, n, nspec):
    #Two mode squeezer
    #for off-diagonal
    m1=mode1*nspec*2+1
    m2=mode2*nspec*2+1
    #for diagonal
    m1s=mode1*nspec
    m2s=mode2*nspec
    c2=Matrix(n,n, lambda i,j: c(xi[i%nspec]) if (i==j)and( 
        ((m1s<=i)and(i<(m1s+nspec)))or((m2s<=i)and(i<(m2s+nspec))) ) else 1 if i==j else 0 )
    s2=Matrix(n,n, lambda i,j: s(xi[i%nspec]) if (((i+j)==m1)or((i+j)==m2))and(abs(i-j)==1) else 0)
    return c2,s2 

def matsubs(mat,aold,anew,bold,bnew,cold,cnew):
    temp1=mat.subs(aold,anew)
    temp2=temp1.subs(bold,bnew)
    temp3=temp2.subs(cold,cnew)
    return temp3

#spatial dim
nspace=4
print('Spatial modes =', nspace)

#spectral dim 
nspectral=2
print('Spectral modes =', nspectral)

#make total dimension
n=nspace*nspectral
print('Total dim =', n)

# 2 mode squeezer on modes 0 & 1
mode1sq,mode2sq= 0,1
print('\n Two mode squeezing on modes, ', mode1sq,',', mode2sq)
# 2 mode sq on modes 2 & 3
mode3sq,mode4sq=2,3
print('Two mode squeezer on modes, ', mode3sq, ',', mode4sq)
# beamsplitter spatial modes
bsmode1, bsmode2= 1,2
print('Beamsplitter on modes, ', bsmode1, ',', bsmode2)

numofsqs =2 

#Show matrices?
showmat=0

#define symbols for modes, 2 spatial, 2 spectral 
a=symbols('a:%d:%d' % (nspace, nspectral)) 
#define symbols for squeezing
xi=symbols('xi:100')
#cosh & sinh placeholders
c,s = symbols('c s')

#Make optical modes
ablk=MatrixSymbol('a',n,1)
#build block matrix for modes
bmodes=BlockMatrix([[ablk],[conjugate(ablk)]])
#make a size n vector for a
modematrix=Matrix(n,1, lambda i,j: a[i] )
#substitution for a modes
modes=matsubs(bmodes,ablk,modematrix,0,0,0,0)
print('\n Mode matrix, a_ij (spatial,spectral)')
pprint(bmodes)

# make squeezing matrices
#Single mode squeezer
c1,s1 = makesinglesq(n, nspectral)
#Two mode squeezer
c12,s12 = maketwomodesq(mode1=mode1sq,mode2=mode2sq, n=n, nspec=nspectral)
c34,s34= maketwomodesq(mode1=mode3sq,mode2=mode4sq, n=n, nspec=nspectral)


#block form for symplectic matrix
al = MatrixSymbol('alpha', n,n)
be = MatrixSymbol('beta',n,n)

#make block form for symplectic matrix
block=BlockMatrix([[al,be],
                [conjugate(be), conjugate(al)]])

print('\n Symplectic matrix')
pprint(block)

ms01=Matrix(matsubs(block,al,c12,be,s12,0,0))
ms23=Matrix(matsubs(block,al,c34,be,s34,0,0))

print('\n Function for debugging')

print('Beamsplitter')
beamspace= Matrix(makebs(bsmode1,bsmode2,nspace,pi/4))
beamsplitter=TensorProduct(beamspace,eye(nspectral))

mbs=Matrix(matsubs(block,al,beamsplitter,be,zeros(n),0,0))

blockbs=block 

#Do mode transformation
transform=block_collapse(blockbs*bmodes)
modeout=matsubs(transform,ablk,modematrix,al,beamsplitter,be,zeros(n))

print('\n Beamsplitter transformation on spatial modes 0 & 1')
if showmat==1:
    pprint(modeout)

print('\n Beamsplitter transforms modes as, ')
pprint(relational.Eq(Matrix(modes[0:n]),Matrix(modeout[0:n])))


#Do mode transformation
transform=block_collapse(block*bmodes)

#Substitute for real this time
print('\n Then Substitute')
finalmatssq= matsubs(transform,ablk,modematrix,al,c1,be,s1)
finalmattmsq=matsubs(transform,ablk,modematrix,al,c12,be,s12)

print('\n Single mode squeezing')
print('\n Print final matrix after subs')
if showmat==1:
    pprint(finalmatssq)

print('\n Single mode squeezing transforms as, ')
pprint(relational.Eq(Matrix(modes[0:n]),Matrix(finalmatssq[0:n])))

print('\nTwo mode squeezing')
print('\n Print final matrix after subs')
if showmat==1:
    pprint(finalmattmsq)

print('\n Two mode squeezing transforms as, ')
outputmodes=Matrix(modes[0:n])
pprint(relational.Eq(Matrix(modes[0:n]),Matrix(finalmattmsq[0:n])))

#############################################################################
############################################################################
#
#                                END OF TESTING 
#
#############################################################################

# beamsplitter spatial modes
print('Beamsplitter on modes, ', bsmode1, ',', bsmode2)
if showmat==1:
    pprint(mbs)

# 2 mode sq on modes 2 & 3
print('Two mode squeezer on modes, ', mode3sq, ',', mode4sq)
if showmat==1:
    pprint(ms23)

# 2 mode squeezer on modes 0 & 1
print('\n Two mode squeezing on modes, ', mode1sq,',', mode2sq)
if showmat==1:
    pprint(ms01)

transform=mbs*ms23*ms01

modetransform=(transform*modes)

pprint(relational.Eq(Matrix(modes[0:n]),Matrix(modetransform[0:n])))
