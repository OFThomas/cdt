from sympy import *
from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum import TensorProduct 
from mpmath import mp
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

def makeps(mode1,phaseangle):
    ##Phase shifter
    m1=mode1*nspectral
    phasespace=Matrix(n,n, lambda i,j: exp(sqrt(-1)*omega[i%nspectral]) if ((i==j)and((m1<=i)and(i<m1+nspectral))) else 1 if(i==j) else 0)
    # symplectic phase shifer
    mps=Matrix(matsubs(block,al,phasespace,be,zeros(n),0,0))
    return mps

def makebs(mode1,mode2,theta):
    #Make beamsplitter 
    beamspace= Matrix(nspace,nspace, lambda i,j: cos(theta) if ((i==j)and((i==mode1)or(i==mode2))) 
            else  sqrt(-1)*sin(theta) if (i+j)==(mode1+mode2)and((i==mode1)or(i==mode2)) 
            else 1 if (i==j)and((i!=mode1)or(i!=mode2)) else 0)
    beamsplitter=TensorProduct(beamspace,eye(nspectral))    
    # symplectic beamsplitter
    mbs=Matrix(matsubs(block,al,beamsplitter,be,zeros(n),0,0))
    return mbs

def makesinglesq(n, nspec):
    #Single mode squeezer
    c1=Matrix(n,n, lambda i,j: c(xi[i%nspec]) if i==j else 0)
    s1=Matrix(n,n, lambda i,j: s(xi[i%nspec]) if i==j else 0)
    return c1, s1
    
def makesq(mode1,mode2, sqparam):
    #Two mode squeezer
    #for off-diagonal
    nspec=nspectral
    m1=mode1*nspec*2+1
    m2=mode2*nspec*2+1
    #for diagonal
    m1s=mode1*nspec
    m2s=mode2*nspec
    #diagonal part cosh on modes, I elsewhere
    c2=Matrix(n,n, lambda i,j: c(sqparam[i]) if (i==j)and( 
        ((m1s<=i)and(i<(m1s+nspec)))or((m2s<=i)and(i<(m2s+nspec))) ) else 1 if i==j else 0 )
    #off diag sinh on modes, zero else
    s2=Matrix(n,n, lambda i,j: s(sqparam[i]) if  (i==n-1-j)and( 
        ((m1s<=i)and(i<(m1s+nspec)))or((m2s<=i)and(i<(m2s+nspec))) )  else 0)
    # symplectic sqs
    ms01=Matrix(matsubs(block,al,c2,be,s2,0,0))
    return ms01

def matsubs(mat,aold,anew,bold,bnew,cold,cnew):
    temp1=mat.subs(aold,anew)
    temp2=temp1.subs(bold,bnew)
    temp3=temp2.subs(cold,cnew)
    return temp3

def justdoitplease(transform, modes):
    modetransform=transform*modes
    pprint(relational.Eq(Matrix(modes[0:n]),Matrix(modetransform[0:n])))
    return 

def takagi_for_unitary(A):
    ### takes a unitary matrix A such that A^T = A, ###
    ###    returns unitary U such that A = U U^T    ###
    N = A.shape[0]
    U = np.eye(N)

    while A.shape[0] > 0:   
        n = A.shape[0]
        # constructing a vector v such that A*ybar = y
        x = np.random.rand(n)+1j*np.random.rand(n)
        x /= np.sqrt(np.dot(x.conj(),x))
        Axbar = np.tensordot(A,x.conj(),(1,0))
        mu = np.dot(x.conj(),Axbar)

        if np.abs(mu)>1-1e-8:
            y = mu**.5*x
        else:
            y = Axbar + x; y /= np.sqrt(np.dot(y.conj(),y))

        # constructing a unitary V with first column = y
        V = np.random.random((n,n))+1j*np.random.random((n,n))
        V[:,0] = y
        V,__ = np.linalg.qr(V)
        # update A
        A = np.tensordot(np.tensordot(V.conj().T,A,(1,0)),V.conj(),(1,0))[1:,1:]
        # update U
        Vi = np.eye(N,dtype=complex)
        Vi[N-n:,N-n:] = V
        U = np.tensordot(U,Vi,(1,0))

    return U

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
xi=[None]*n
for i in range(0,n):
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
msq1=makesq(sq1_mode1,sq1_mode2, xi)
pprint(msq1)

#make tmsq on modes 2,3
msq2=makesq(sq2_mode1,sq2_mode2, xi)
pprint(msq2)

#make phase shifter
mps=makeps(phasemode,phaseangle)
pprint(mps)

#make beamsplitter 
mbs=makebs(bsmode1,bsmode2,bsangle)
pprint(mbs)

#do mode transformation
justdoitplease(mbs*mps*msq2*msq1,modes)

print('Worked?')

pprint(msq1)
S1=msq1[range(sq1_mode1,sq1_mode2+1),range(n+sq1_mode1,n+sq1_mode2+1)]
pprint(S1)

unitary=diag(p, eye(2))
#pprint(unitary)
u1=TensorProduct(eye(2),unitary)
#pprint(u1)

#pprint(msq1)
#pprint(u1.inv()*msq1*u1)

#p,q =symbols('p q')
#p=2+sqrt(-1)*3
#q=3-sqrt(-1)

print('\nTesting diag structure on symplectic matrices\n')
#sym=Matrix([[p,q],[conjugate(q),conjugate(p)]])
#pprint(sym)

print('\nnumerical Svd')
#u,s,v = mp.svd(sym)
#pprint(u)
#pprint(s)
#pprint(v)

print()
print('Columns of U \neigenvects of sym*conj(sym)')
#pprint((sym*conjugate(sym)).eigenvects())

print('Columns of V\neigenvects of conj(sym)*sym')
#pprint((conjugate(sym)*sym).eigenvects())


"""
######################## Substitution for numerical values
xi[0]=0
xi[1]='blue'

#update matrices and mode transformation

   
#make tmsq on modes 0,1 squeezing xi
msq1=makesq(sq1_mode1,sq1_mode2, xi)
#make tmsq on modes 2,3
msq2=makesq(sq2_mode1,sq2_mode2, xi)
#make phase shifter
mps=makeps(phasemode,phaseangle)
#make beamsplitter 
mbs=makebs(bsmode1,bsmode2,bsangle)

#do mode transformation
justdoitplease(mbs*mps*msq2*msq1,modes)

"""


#end

