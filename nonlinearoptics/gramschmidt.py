from sympy import *

l=[Matrix([1,0,0]), Matrix([1/3,1/3,1/3])]
out=GramSchmidt(l,True)
pprint(out)
