from sympy import *
from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum import TensorProduct 
from sympy.abc import alpha, beta, xi, zeta
init_printing(use_unicode=True)

def makemodes():
    #Make optical modes
    ablk=MatrixSymbol('a',n,1)
    #build block matrix for modes
    bmodes=BlockMatrix([[ablk],[conjugate(ablk)]])
    #make a size n vector for a
    modematrix=Matrix(n,1, lambda i,j: a[i] )
    #substitution for a modes
    modes=matsubs(bmodes,ablk,modematrix,0,0,0,0)
    return modes

def makeps(mode1,n,nspec,theta):
    #Phase shifter
    m1=mode1*nspec
    passive=Matrix(n,n, lambda i,j: exp(sqrt(-1)*omega[i%nspec]) if ((i==j)and((m1<=i)and(i<m1+nspec))) else 1 if(i==j) else 0)
    return passive

def makemps(m1,phaseangle):
    #Make phase shifter 
    phasespace=Matrix(makeps(m1,n,nspectral,phaseangle))
    # symplectic phase shifer
    mps=Matrix(matsubs(block,al,phasespace,be,zeros(n),0,0))
    return mps

def makebs(mode1,mode2, n, theta):
    #beamsplitter
    passive=Matrix(n,n, lambda i,j: cos(theta) if ((i==j)and((i==mode1)or(i==mode2))) 
            else  sqrt(-1)*sin(theta) if (i+j)==(mode1+mode2)and((i==mode1)or(i==mode2)) 
            else 1 if (i==j)and((i!=mode1)or(i!=mode2)) else 0)
    return passive

def makembs(m1,m2,angle):
    #Make beamsplitter 
    beamspace= Matrix(makebs(m1,m2,nspace,angle))
    beamsplitter=TensorProduct(beamspace,eye(nspectral))    
    # symplectic beamsplitter
    mbs=Matrix(matsubs(block,al,beamsplitter,be,zeros(n),0,0))
    return mbs

def makesinglesq(n, nspec):
    #Single mode squeezer
    c1=Matrix(n,n, lambda i,j: c(xi[i%nspec]) if i==j else 0)
    s1=Matrix(n,n, lambda i,j: s(xi[i%nspec]) if i==j else 0)
    return c1, s1
    
def maketwomodesq(mode1,mode2, n, nspec, sqparam):
    #Two mode squeezer
    #for off-diagonal
    m1=mode1*nspec*2+1
    m2=mode2*nspec*2+1
    #for diagonal
    m1s=mode1*nspec
    m2s=mode2*nspec
    #diagonal part cosh on modes, I elsewhere
    c2=Matrix(n,n, lambda i,j: c(sqparam[i%nspec]) if (i==j)and( 
        ((m1s<=i)and(i<(m1s+nspec)))or((m2s<=i)and(i<(m2s+nspec))) ) else 1 if i==j else 0 )
    #off diag sinh on modes, zero else
    s2=Matrix(n,n, lambda i,j: s(sqparam[i%nspec]) if (((i+j)==m1)or((i+j)==m2))and(abs(i-j)==1) else 0)
    return c2,s2 

def makemsq(m1,m2, xi):
    #Two mode squeezer
    c12,s12 = maketwomodesq(mode1=m1,mode2=m2, n=n, nspec=nspectral, sqparam=xi)
    # symplectic sqs
    ms01=Matrix(matsubs(block,al,c12,be,s12,0,0))
    return ms01

def matsubs(mat,aold,anew,bold,bnew,cold,cnew):
    temp1=mat.subs(aold,anew)
    temp2=temp1.subs(bold,bnew)
    temp3=temp2.subs(cold,cnew)
    return temp3

def doitplease():

    #Make optical modes
    ablk=MatrixSymbol('a',n,1)
    #build block matrix for modes
    bmodes=BlockMatrix([[ablk],[conjugate(ablk)]])
    #make a size n vector for a
    modematrix=Matrix(n,1, lambda i,j: a[i] )
    #substitution for a modes
    modes=matsubs(bmodes,ablk,modematrix,0,0,0,0)

    #block form for symplectic matrix
    al = MatrixSymbol('alpha', n,n)
    be = MatrixSymbol('beta',n,n)

    #make block form for symplectic matrix
    block=BlockMatrix([[al,be],
                 [conjugate(be), conjugate(al)]])

    # make squeezing matrices
    #Single mode squeezer
    c1,s1 = makesinglesq(n, nspectral)
    #Two mode squeezer
    c12,s12 = maketwomodesq(mode1=sq1_mode1,mode2=sq1_mode2, n=n, nspec=nspectral, sqparam=xi)
    c34,s34= maketwomodesq(mode1=sq2_mode1,mode2=sq2_mode2, n=n, nspec=nspectral, sqparam=xi)
 
    # symplectic sqs
    ms01=Matrix(matsubs(block,al,c12,be,s12,0,0))
    ms23=Matrix(matsubs(block,al,c34,be,s34,0,0))

    #Make beamsplitter 
    beamspace= Matrix(makebs(bsmode1,bsmode2,nspace,pi/4))
    beamsplitter=TensorProduct(beamspace,eye(nspectral))
    
    # symplectic beamsplitter
    mbs=Matrix(matsubs(block,al,beamsplitter,be,zeros(n),0,0))

    #Make phase shifter 
    phasespace=Matrix(makeps(phasemode,n,nspectral,pi/2))
    
    # symplectic phase shifer
    mps=Matrix(matsubs(block,al,phasespace,be,zeros(n),0,0))
    # return symplectic transform and modes 
    print('\nms01, ms23, mps, mbs\n')
    return mbs*mps*ms23*ms01, modes


def justdoitplease(transform, modes):
    modetransform=transform*modes
    pprint(relational.Eq(Matrix(modes[0:n]),Matrix(modetransform[0:n])))
    return 

################################## start of program

#spatial dim
nspace=4

#spectral dim 
nspectral=2

#make total dimension
n=nspace*nspectral

################# Specify modes for unitaries
# 2 mode squeezer on modes 0 & 1
sq1_mode1,sq1_mode2= 0,1

# 2 mode sq on modes 2 & 3
sq2_mode1,sq2_mode2=2,3

# Phase shifter
phasemode=2
phaseangle=pi/2

# beamsplitter spatial modes
bsmode1, bsmode2= 1,2
bsangle=pi/4

# phase shift
omega=symbols('omega:100')

print('Spatial modes =', nspace)
print('Spectral modes =', nspectral)
print('Total dim =', n)
print('\n Two mode squeezing on modes, ', sq1_mode1,',', sq1_mode2)
print('Two mode squeezer on modes, ', sq2_mode1, ',', sq2_mode2)
print('Phase shift on modes. ', phasemode, 'phase angle', phaseangle)
print('Beamsplitter on modes, ', bsmode1, ',', bsmode2, 'angle= ', bsangle)

########### Define squeezing symbols
xi=[None]*2
for i in range(0,2):
    xi[i]=symbols('xi%d' % (i))

#cosh & sinh placeholders
c,s = symbols('c s')

########### Define optical modes
a=[None]*n
for i in range(0,nspace):
    for j in range(0,nspectral):
        a[i*(nspectral)+j]=symbols('a%d%d' % (i, j))

################## Make stuff happen!
modes=makemodes()

#block form for symplectic matrix
al = MatrixSymbol('alpha', n,n)
be = MatrixSymbol('beta',n,n)

#make block form for symplectic matrix
block=BlockMatrix([[al,be],
    [conjugate(be), conjugate(al)]])
   
#make tmsq on modes 0,1 squeezing xi
msq1=makemsq(sq1_mode1,sq1_mode2, xi)
#make tmsq on modes 2,3
msq2=makemsq(sq2_mode1,sq2_mode2, xi)
#make phase shifter
mps=makemps(phasemode,phaseangle)
#make beamsplitter 
mbs=makembs(bsmode1,bsmode2,bsangle)

#do mode transformation
justdoitplease(mbs*mps*msq2*msq1,modes)

print('Worked?')

######################## Substitution for numerical values
xi[0]=0
xi[1]='blue'

#update matrices and mode transformation

   
#make tmsq on modes 0,1 squeezing xi
msq1=makemsq(sq1_mode1,sq1_mode2, xi)
#make tmsq on modes 2,3
msq2=makemsq(sq2_mode1,sq2_mode2, xi)
#make phase shifter
mps=makemps(phasemode,phaseangle)
#make beamsplitter 
mbs=makembs(bsmode1,bsmode2,bsangle)

#do mode transformation
justdoitplease(mbs*mps*msq2*msq1,modes)




#end

