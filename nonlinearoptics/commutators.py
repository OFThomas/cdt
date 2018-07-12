from sympy.physics.secondquant import Commutator as Com
from sympy import *
from sympy.physics.secondquant import B, Dagger  

a=[None]*2

for i in range(0,len(a)):
    a[i]= B('%d' % (i))


pprint(Com(Dagger(a[1]),(a[1])))

print((a[0]))
# pprint(com.doit())
