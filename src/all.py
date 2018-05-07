from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

excase=int(raw_input('Enter the example number \n'))

###########################################################################################################

def ex1():
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
    return 1
############################################################################################################

def ex2():
    # -*- coding: utf-8 -*-
    """
    Created on Fri Dec 16 09:30:30 2011
    Python Batman Equation
    @author: Trae Blain
    """
    from numpy import sqrt
    from numpy import meshgrid
    from numpy import arange

    xs = arange(-7.25, 7.25, 0.01)
    ys = arange(-5, 5, 0.01)
    x, y = meshgrid(xs, ys)

    eq1 = ((x/7)**2*sqrt(abs(abs(x)-3)/(abs(x)-3))+(y/3)**2*sqrt(abs(y+3/7*sqrt(33))/(y+3/7*sqrt(33)))-1)
    eq2 = (abs(x/2)-((3*sqrt(33)-7)/112)*x**2-3+sqrt(1-(abs(abs(x)-2)-1)**2)-y)
    eq3 = (9*sqrt(abs((abs(x)-1)*(abs(x)-.75))/((1-abs(x))*(abs(x)-.75)))-8*abs(x)-y)
    eq4 = (3*abs(x)+.75*sqrt(abs((abs(x)-.75)*(abs(x)-.5))/((.75-abs(x))*(abs(x)-.5)))-y)
    eq5 = (2.25*sqrt(abs((x-.5)*(x+.5))/((.5-x)*(.5+x)))-y)
    eq6 = (6*sqrt(10)/7+(1.5-.5*abs(x))*sqrt(abs(abs(x)-1)/(abs(x)-1))-(6*sqrt(10)/14)*sqrt(4-(abs(x)-1)**2)-y)

    #eq1 = ((x/7.0)**2.0*sqrt(abs(abs(x)-3.0)/(abs(x)-3.0))+(y/3.0)**2.0*sqrt(abs(y+3.0/7.0*sqrt(33.0))/(y+3.0/7.0*sqrt(33.0)))-1.0)

    for f in [eq1,eq2,eq3,eq4,eq5,eq6]:
        plt.contour(x, y, f, [0])

    plt.title("I'M BATMAN.")

    #plt.savefig('ex2.png', format='png')
    plt.show()
    return 1
###########################################################################################################

def ex3():
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
    return 1
############################################################################

def ex4():
    import numpy as np
    import matplotlib.mlab as mlab
    import matplotlib.pyplot as plt

    mu, sigma = 100, 15
    x = mu + sigma*np.random.randn(10000)

    # the histogram of the data
    n, bins, patches = plt.hist(x, 50, normed=1, facecolor='green', alpha=0.75)

    # add a 'best fit' line
    y = mlab.normpdf( bins, mu, sigma)
    l = plt.plot(bins, y, 'r--', linewidth=1)
    
    plt.xlabel('Smarts')
    plt.ylabel('Probability')
    plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')

    plt.axis([40, 160, 0, 0.03])
    plt.grid(True)

    plt.show()
    #plt.savefig('ex4.png', format='png')
    return 1
#############################################################################

def ex5():
    """
    =================
    Multiple subplots
    =================
    Simple demo with multiple subplots.
    """
    import numpy as np
    import matplotlib.pyplot as plt

    x1 = np.linspace(0.0, 5.0)
    x2 = np.linspace(0.0, 2.0)

    y1 = np.cos(2 * np.pi * x1) * np.exp(-x1)
    y2 = np.cos(2 * np.pi * x2)

    plt.subplot(2, 1, 1)
    plt.plot(x1, y1, 'o-')
    plt.title('A tale of 2 subplots')
    plt.ylabel('Damped oscillation')

    plt.subplot(2, 1, 2)
    plt.plot(x2, y2, '.-')
    plt.xlabel('time (s)')
    plt.ylabel('Undamped')

    plt.show()
    #plt.savefig('ex5.png', format='png')
    return 1
###############################################################################

def ex6():
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, PathPatch
    from matplotlib.path import Path
    from matplotlib.transforms import Affine2D
    import numpy as np

    # Fixing random state for reproducibility
    np.random.seed(19680801)

    r = np.random.rand(50)
    t = np.random.rand(50) * np.pi * 2.0
    x = r * np.cos(t)
    y = r * np.sin(t)

    fig, ax = plt.subplots(figsize=(6, 6))
    circle = Circle((0, 0), 1, facecolor='none',
                edgecolor=(0, 0.8, 0.8), linewidth=3, alpha=0.5)
    ax.add_patch(circle)

    im = plt.imshow(np.random.random((100, 100)),
                    origin='lower', cmap=cm.winter,
                    interpolation='spline36',
                    extent=([-1, 1, -1, 1]))
    im.set_clip_path(circle)

    plt.plot(x, y, 'o', color=(0.9, 0.9, 1.0), alpha=0.8)

    # Dolphin from OpenClipart library by Andy Fitzsimon
    #       <cc:License rdf:about="http://web.resource.org/cc/PublicDomain">
    #         <cc:permits rdf:resource="http://web.resource.org/cc/Reproduction"/>
    #         <cc:permits rdf:resource="http://web.resource.org/cc/Distribution"/>
    #         <cc:permits rdf:resource="http://web.resource.org/cc/DerivativeWorks"/>
    #       </cc:License>

    dolphin = """
    M -0.59739425,160.18173 C -0.62740401,160.18885 -0.57867129,160.11183
    -0.57867129,160.11183 C -0.57867129,160.11183 -0.5438361,159.89315
    -0.39514638,159.81496 C -0.24645668,159.73678 -0.18316813,159.71981
    -0.18316813,159.71981 C -0.18316813,159.71981 -0.10322971,159.58124
    -0.057804323,159.58725 C -0.029723983,159.58913 -0.061841603,159.60356
    -0.071265813,159.62815 C -0.080250183,159.65325 -0.082918513,159.70554
    -0.061841203,159.71248 C -0.040763903,159.7194 -0.0066711426,159.71091
    0.077336307,159.73612 C 0.16879567,159.76377 0.28380306,159.86448
    0.31516668,159.91533 C 0.3465303,159.96618 0.5011127,160.1771
    0.5011127,160.1771 C 0.63668998,160.19238 0.67763022,160.31259
    0.66556395,160.32668 C 0.65339985,160.34212 0.66350443,160.33642
    0.64907098,160.33088 C 0.63463742,160.32533 0.61309688,160.297
    0.5789627,160.29339 C 0.54348657,160.28968 0.52329693,160.27674
    0.50728856,160.27737 C 0.49060916,160.27795 0.48965803,160.31565
    0.46114204,160.33673 C 0.43329696,160.35786 0.4570711,160.39871
    0.43309565,160.40685 C 0.4105108,160.41442 0.39416631,160.33027
    0.3954995,160.2935 C 0.39683269,160.25672 0.43807996,160.21522
    0.44567915,160.19734 C 0.45327833,160.17946 0.27946869,159.9424
    -0.061852613,159.99845 C -0.083965233,160.0427 -0.26176109,160.06683
    -0.26176109,160.06683 C -0.30127962,160.07028 -0.21167141,160.09731
    -0.24649368,160.1011 C -0.32642366,160.11569 -0.34521187,160.06895
    -0.40622293,160.0819 C -0.467234,160.09485 -0.56738444,160.17461
    -0.59739425,160.18173
    """

    vertices = []
    codes = []
    parts = dolphin.split()
    i = 0
    code_map = {
        'M': (Path.MOVETO, 1),
    	'C': (Path.CURVE4, 3),
    	'L': (Path.LINETO, 1)}

    while i < len(parts):
        code = parts[i]
    	path_code, npoints = code_map[code]
    	codes.extend([path_code] * npoints)
    	vertices.extend([[float(x) for x in y.split(',')] for y in
                parts[i + 1:i + npoints + 1]])
    	i += npoints + 1

    vertices = np.array(vertices, float)
    vertices[:, 1] -= 160

    dolphin_path = Path(vertices, codes)
    dolphin_patch = PathPatch(dolphin_path, facecolor=(0.6, 0.6, 0.6),
                             edgecolor=(0.0, 0.0, 0.0))

    ax.add_patch(dolphin_patch)

    vertices = Affine2D().rotate_deg(60).transform(vertices)
    dolphin_path2 = Path(vertices, codes)
    dolphin_patch2 = PathPatch(dolphin_path2, facecolor=(0.5, 0.5, 0.5),
                            edgecolor=(0.0, 0.0, 0.0))                    
    ax.add_patch(dolphin_patch2)

    #plt.savefig('ex6.png', format='png')
    plt.show()
    return 1
##################################################################################

if(excase==1):
    ex1()
elif(excase==2):
    ex2()
elif(excase==3):
    ex3()
elif(excase==4):
    ex4()
elif(excase==5):
    ex5()
elif(excase==6):
    ex6()
else:
    print 'Not an int between 1 & 6'


