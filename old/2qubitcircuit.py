import math
from qutip import *

n=4
qc = QubitCircuit(n)
qc.add_gate("X", 2,[0,1])
#qc.png

ulist=qc.propagators()
u=gate_sequence_product(ulist)
print(u)


