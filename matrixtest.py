from sympy import *
init_printing(use_unicode=True)

a = symbols('a')

x=Matrix([[0,1],[1,0]])
pprint (x)

print x.eigenvals()
pprint(x.eigenvects())

p, d = x.diagonalize()
pprint(p)
pprint(d)

pprint(p*d*p**-1)





