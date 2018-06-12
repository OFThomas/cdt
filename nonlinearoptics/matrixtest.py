from sympy import *
init_printing(use_unicode=True)

c1=MatrixSymbol('c1',2,2)
s1=MatrixSymbol('s1',2,2)

#cbig=MatrixSymbol('c', 2,2)
#sbig=MatrixSymbol('s', 2,2)
c,s = symbols('c,s')
cbig, sbig = symbols('c s')


a,b,c,d = symbols('a,b,c,d')
e,f,g,h = symbols('e,f,g,h')

#Squeezing mode size
sqsize=2

#total dim
n=2

#make symbolic matrix
smsq=Matrix(n,n, lambda i,j: cbig if i==j else sbig if i+j==n-1 else 0)
print('smsq')
pprint(smsq)

#single mode squeezer
c=Matrix(sqsize,sqsize, lambda i,j: a )
s=Matrix(sqsize,sqsize, lambda i,j: b)

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


"""
Eq(a,1)

pprint (x)

print(x.eigenvals())
pprint(x.eigenvects())

p, d = x.diagonalize()
pprint(p)
pprint(d)

pprint(p*d*p**-1)

"""


