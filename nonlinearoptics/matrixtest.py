from sympy import *
from sympy.physics.quantum.dagger import Dagger
from sympy.abc import alpha, beta, xi, zeta
init_printing(use_unicode=True)

#single mode squeezer 
def makesinglesq(n):
    #Single mode squeezer
    c1=Matrix(n,n, lambda i,j: c(xi[i%2]) if i==j else 0)
    s1=Matrix(n,n, lambda i,j: s(xi[i%2]) if i==j else 0)
    return c1, s1
    
def maketwomodesq(n):
    #Two mode squeezer
    c2=Matrix(n,n, lambda i,j: c(xi[i%2]) if i==j else 0 )
    s2=Matrix(n,n, lambda i,j: s(xi[i%2]) if (i+j)%n==1 else 0)
    return c2,s2 

def matsubs(mat,aold,anew,bold,bnew,cold,cnew):
    temp1=mat.subs(aold,anew)
    temp2=temp1.subs(bold,bnew)
    temp3=temp2.subs(cold,cnew)
    return temp3

#spatial dim
nspace=2
print('Spatial modes =', nspace)

#spectral dim 
nspectral=2
print('Spectral modes =', nspectral)

#make total dimension
n=nspace*nspectral
print('Total dim =', n)

#define symbols for modes, 2 spatial, 2 spectral 
a=symbols('a:2:2')
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

#Single mode squeezer
c1,s1 = makesinglesq(n)

#Two mode squeezer
c2,s2 = maketwomodesq(n)

#block form for symplectic matrix
al = MatrixSymbol('alpha', n,n)
be = MatrixSymbol('beta',n,n)

#make block form for symplectic matrix
block=BlockMatrix([[al,be],
                [conjugate(be), conjugate(al)]])

print('\n Symplectic matrix')
pprint(block)

print('\n Function for debugging')

#Do mode transformation
transform=block_collapse(block*bmodes)

#Substitute for real this time
print('\n Then Substitute')
finalmatssq= matsubs(transform,ablk,modematrix,al,c1,be,s1)
finalmattmsq=matsubs(transform,ablk,modematrix,al,c2,be,s2)

print('\n Single mode squeezing')
print('\n Print final matrix after subs')
pprint(finalmatssq)

print('\n Which transforms the modes as, ')
pprint(relational.Eq(Matrix(modes[0:n]),Matrix(finalmatssq[0:n])))

print('\nTwo mode squeezing')
print('\n Print final matrix after subs')
pprint(finalmattmsq)

print('\n Two mode squeezing transforms as, ')
pprint(relational.Eq(Matrix(modes[0:n]),Matrix(finalmattmsq[0:n])))


