from sympy import *
from sympy.abc import alpha, beta
init_printing(use_unicode=True)


def makesinglesq(size, sqcoeff):
    for i in range(1,size):
        return 1

c1=MatrixSymbol('c1',2,2)
s1=MatrixSymbol('s1',2,2)

#cbig=MatrixSymbol('c', 2,2)
#sbig=MatrixSymbol('s', 2,2)
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

#for i in range(1,2):
 #   al[i,i]=1
mata=Matrix([[1,0],[0,1]])
matb=Matrix([[0,1],[1,0]])

block=BlockMatrix([[al,be],
                [conjugate(be), conjugate(al)]])
pprint(block)
print()
pprint(Matrix(block))

temp=block.subs(al,mata)
end=Matrix(temp.subs(be,matb))
#end=ImmutableMatrix(tend)
pprint(end)
print(type(end))
#p,d=end.diagonalize()
#mend=end.clone()
mend=end.evalf()
print(type(mend))
#pprint(mend.eigenvectors())
#pprint(block_collapse(block.diagonalize()))
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


