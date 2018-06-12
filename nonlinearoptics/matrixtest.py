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

mata=Matrix([[1,0],[0,1]])
matb=Matrix([[0,1],[1,0]])

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
end=Matrix(temp.subs(be,matb))
pprint(end)
print(type(end))
#p,d=end.diagonalize()
#mend=end.clone()
mend=end.evalf()
print(type(mend))
print(type(mata))
pprint(end)

print('\nEigenvectors\n')
pprint(end.eigenvects())

p,d = end.diagonalize()

pprint(p)
print('\n diag matrix using PDP^-1\n')
pprint(d)

pprint(exp(end))


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


