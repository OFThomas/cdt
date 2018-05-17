import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def semicirc(x):
    return np.sqrt(1-x**2)

def ellipse(x):
    return 3*semicirc(x/7)

def shoulders(x):
    return 4.2 - 0.5*x - 2.8*semicirc(0.5*x-0.5)

def bot(x):
    return semicirc(abs(2-x)-1) - x**2/11 + 0.5*x -3

def xl(min,max):
    return np.linspace(min,max)
#sf(x) = sqrt(1-x^2)                           # semicircle
#ef(x) = 3*sf(x/7)                             # ellipse
#sh(x) = 4.2 - .5*x -2.8*sf(.5*x -.5)          # shoulders
#bf(x) = sf(abs(2 - x) - 1) - x^2/11 + .5*x -3 # bottom

cl    = [(0,1.7), (.5,1.7), (.8,2.6), (1,.9)] # cowl right
cl2   = [(0,1.7), (-.5,1.7), (-.8,2.6), (-1,.9)]               # cowl left

def p(f,xmin,xmax):
#    "symmetric plot across y-axis"
    p1 = plot(f,xmin,xmax)
    p2 = plot(lambda x:f(-x),-xmax,-xmin)
    return p1 + p2

#p(ellipse,3,7) + p(-ellipse,4,7) + p(shoulders,1,3) + p(bot,0,4) + line(cl) +line(cl2)
x = np.linspace(-10,10)

y=ellipse(xl(3,7)) -ellipse(xl(4,7)) + shoulders(xl(1,3)) + bot(xl(0,4)) 

plt.plot(x, y) 

plt.show()

