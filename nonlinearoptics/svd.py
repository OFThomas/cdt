from sympy import *
init_printing(use_unicode=True)

c,s = symbols('c s')

n=4

# lambda c on diagonals, s on anti diag
diagc=Matrix(n,n, lambda i,j: c if i==j else s if i+j==n-1 else 0)

pprint(diagc)



