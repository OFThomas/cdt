import numpy as np
import matplotlib.pyplot as plt

# evenly sampled points 
t = np.linspace(0, 5)

# red dashes, blue squares and green triangles
plt.plot(t, t, 'r--', t, t**2, 'bs', t, t**3, 'g^')

#labels
plt.ylabel('Y data /units')
plt.xlabel('X data /units')

plt.title('3 functions on one graph!')

#plt.savefig('ex1.png', format='png')
plt.show()
