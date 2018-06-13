from sympy import *
from sympy.physics.quantum.dagger import Dagger
from sympy.abc import alpha, beta, xi, zeta
init_printing(use_unicode=True)

def makesinglesq(size, sqcoeff):
    for i in range(1,size):
        return 1

xi1, xi2 = symbols('xi_1 xi_2')

pprint(xi1)

c1=MatrixSymbol('c1',2,2)
s1=MatrixSymbol('s1',2,2)

c,s = symbols('c,s')
cbig, sbig = symbols('c s')

a,b,c,d = symbols('a,b,c,d')
e,f,g,h = symbols('e,f,g,h')

sqcoeffs=[s for i in range(2) ]
print(sqcoeffs)

sqmatrix=makesinglesq(2,sqcoeffs)
pprint(sqmatrix)
#Squeezing mode size
sqsize=2

#total dim
n=2

al = MatrixSymbol('alpha', 2,2)
be = MatrixSymbol('beta',2,2)

mata=Matrix([[xi1,0],[0,xi1]])
matb=Matrix([[0,xi2],[xi2,0]])

#make block form for matrix
block=BlockMatrix([[al,be],
                [conjugate(be), conjugate(al)]])

print('Symplectic matrix\n')
pprint(block)
print()
print('Break the blocks up\n')
mblock=Matrix(block)
pprint(mblock)
print()

#substitute alpha
temp=block.subs(al,mata)
#subs beta
blockend=temp.subs(be,matb)
end=Matrix(blockend)

pprint(end)
print(type(end))
#p,d=end.diagonalize()
#mend=end.clone()
mend=end.evalf()
print(type(mend))
print(type(mata))
pprint(end)

print('\nEigenvectors\n')
pprint(end.eigenvals())

p,d = end.diagonalize()

pprint(p)
print('\n diag matrix using PDP^-1\n')
pprint(d)
print('\n simplify\n')
#pprint(exp(end))


sqmat=Matrix([[cbig,sbig],
            [conjugate(sbig), conjugate(cbig)]])

print('make block form\n')
pprint(sqmat)

sqsize=2
#single mode squeezer
c=Matrix(sqsize,sqsize, lambda i,j: xi1 if i==j else 0 )
s=Matrix(sqsize,sqsize, lambda i,j: xi2 if i+j==sqsize-1 else 0)

test=sqmat.subs({cbig:c, sbig:s})
print('\n subs in blocks\n')
pprint(test)
print(type(test))

print('\n collapse to matrix\n')
testmat=Matrix(test)
pprint(testmat)
print(type(testmat))

pprint(test[0,1])
p,d = test[0,1].diagonalize()
pprint(d)

a=MatrixSymbol('a',2,1)
a1,a2 = symbols('a1 a2')
aamatrix=Matrix([[a1],[a2]])
pprint(aamatrix)

bmodes=BlockMatrix([[a],[conjugate(a)]])
pprint(bmodes)
print('subs in a')
bmodesexp=bmodes.subs(a,aamatrix)
modes=Matrix(bmodesexp)

pprint(modes)

pprint(block_collapse(block*bmodes))
"""
block1=BlockMatrix([[c,s],
                    [conjugate(s),conjugate(c)]])

pprint(block1)
print(type(block1))
pprint(Matrix(block1))
print(type(Matrix(block1)))

upper=block1[0,0]
pprint(upper)
"""

"""
#make symbolic matrix
smsq=Matrix(n,n, lambda i,j: cbig if i==j else sbig if i+j==n-1 else 0)
print('smsq')
pprint(smsq)

#single mode squeezer
c=Matrix(sqsize,sqsize, lambda i,j: 0 )
s=Matrix(sqsize,sqsize, lambda i,j: 1)

testsmsq=Matrix(n,n, lambda i,j: c if i==j else s if i+j==n-1 else 0)
print('testsmsq')
pprint(testsmsq)

print('subs in matrices')
test=smsq.subs({cbig:c, sbig:s})
pprint(test)
smsq=test
print('Print smsq')
pprint(smsq)

#print('eigenvlues, multiplicities')
#print(smsq.eigenvals())
#print('eigenvectors')
#pprint(x.eigenvects())

p, d = smsq.diagonalize()
pprint(p)
pprint(d)

pprint(p*d*p**-1)



Eq(a,1)

pprint (x)

print(x.eigenvals())
pprint(x.eigenvects())

p, d = x.diagonalize()
pprint(p)
pprint(d)

pprint(p*d*p**-1)

"""


