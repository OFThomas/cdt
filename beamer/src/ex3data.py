import numpy as np
import matplotlib.pyplot as plt

def fixy(x,y,coeff, n):
    yfixed = y-(coeff[0]*x + coeff[1])
    return yfixed

#Data file to read in
datafile='data.txt'
datafile=raw_input('enter data file name: \n')
x, y = np.genfromtxt(datafile, unpack=True)
   
#order of fit
n=int(raw_input('enter degrees of freedom for fit: \n'))

#Calculate best fit order n 
#p[0] = ax^2 ,p[1] = bx, p[2] = c
p = np.polyfit(x,y,n)
line=np.polyfit(x,y,1)
    
#don't worry about this line
line[1]=p[2]

#rewrite fit as a function
fit =np.poly1d(p)
#generate smooth points to plot the fit
xp=np.linspace(min(x),max(x), 100)

#PLOT ENERGY DATA
plt.plot(x, y, 'ro', xp, fit(xp), 'r-')

########## Then change title and axis respectively #################

#save graph 
d_name = datafile + '.png'
plt.savefig('ex3'+d_name, format='png')
    
#plot fixed energy of formation
plt.ylabel('Energy, eV')
plt.xlabel('Percentage composition of Mg/Ca')
plt.title('Formation Energy for Composition of Mg/Ca')
    
#save graph
plt.show()
