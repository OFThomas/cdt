from sympy import *
from sympy.physics.quantum.dagger import Dagger
from sympy.abc import alpha, beta, xi, zeta
init_printing(use_unicode=True)

#single mode squeezer 
def makesinglesq(size, sqcoeff):
    for i in range(1,size):
        return 1

def matsubs(mat,aold,anew,bold,bnew,cold,cnew):
    temp1=mat.subs(aold,anew)
    temp2=temp1.subs(bold,bnew)
    temp3=temp2.subs(cold,cnew)
    return temp3

#Squeezing mode size
sqsize=2

#total dim
nspace=2
print('Spatial modes =', nspace)
#spectral dim 
nspectral=2
print('Spectral modes =', nspectral)
#make total dimension
n=nspace*nspectral
print('Total dim =', n)

#define symbols for modes
a=symbols('a:100')

#define symbols for squeezing
xi=symbols('xi:100')

#block form for symplectic matrix
al = MatrixSymbol('alpha', n,n)
be = MatrixSymbol('beta',n,n)

#Make optical modes
ablk=MatrixSymbol('a',n,1)

modematrix=Matrix(n,1, lambda i,j: a[i] )

print('\n Mode matrix')
bmodes=BlockMatrix([[ablk],[conjugate(ablk)]])
modes=matsubs(bmodes,ablk,modematrix,0,0,0,0)
pprint(bmodes)

#single mode squeezer
c=Matrix(n,n, lambda i,j: xi[0] if i==j else 0 )
s=Matrix(n,n, lambda i,j: xi[1] if i+j==n-1 else 0)

#make block form for matrix
block=BlockMatrix([[al,be],
                [conjugate(be), conjugate(al)]])

print('\n Symplectic matrix')
pprint(block)

print('\n Function for debugging')
#Substitute matrix elements in!
#blk=matsubs(block,al,c,be,s,0,0)
#pprint(Matrix(blk))

#Do mode transformation
transform=block_collapse(block*bmodes)

#Substitute for real this time
print('\n Then Substitute')
finalmat = matsubs(transform,ablk,modematrix,al,c,be,s)

print('\n Print final matrix after subs')
pprint(finalmat)
print('\n Which for the input modes evaluates to')
pprint(Matrix(finalmat[0:4]))

print('\n Which transforms the modes as, ')
pprint(relational.Eq(Matrix(modes[0:n]),Matrix(finalmat[0:n])))



