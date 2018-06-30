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
    phasespace=Matrix(n,n, lambda i,j: exp(sqrt(-1)*omega[i%nspectral])
            if ((i==j)and((m1<=i)and(i<m1+nspectral))) else 1 if(i==j) else 0)
    # symplectic phase shifer
    mps=Matrix(matsubs(block,al,phasespace,be,zeros(n),0,0))
    return mps

def makebs(mode1,mode2,theta):
    #Make beamsplitter 
    beamspace= Matrix(nspace,nspace, lambda i,j: cos(theta) 
            if ((i==j)and((i==mode1)or(i==mode2))) 
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
    #for diagonal
    m1s=mode1*nspec
    m2s=mode2*nspec
    #diagonal part cosh on modes, I elsewhere
    c2=Matrix(n,n, lambda i,j: cosh(sqparam[i+m2s]) 
            if (i==j)and((m1s<=i<m1s+nspec)or(m2s<=i<m2s+nspec))
            else 1 if i==j else 0 )
    #off diag sinh on modes, zero else
    s2=Matrix(n,n, lambda i,j: sinh(sqparam[2*i+j-m1s]) 
            if ((m1s<=i<m1s+nspec)or(m2s<=i<m2s+nspec))and
            ((m1s<=j<m1s+nspec)or(m2s<=j<m2s+nspec))
            else 0)
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

def fixfloaterr(matrix):
    imax,jmax=matrix.shape
    newmatrix=matrix
    for j in range(0,jmax):
        for i in range(0,imax):
            if (abs(newmatrix[i,j])<=10**-15):
                newmatrix[i,j]=0
    
    return Matrix(newmatrix)

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
nspectral=1

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
xi=[None]*(n**2)
for i in range(0,n**2):
    xi[i]=symbols('xi%d' % (i))

#modes 0 & 1
xi[1]=1
xi[2]=2

#modes 2 & 3
xi[5]=1
xi[6]=2

#make active s block anti-diag
xi[0]=0
xi[3]=0
xi[4]=0
xi[7]=0

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
S1=msq1[range(sq1_mode1,sq1_mode2+1),range(sq1_mode1+n,sq1_mode2+n+1)]

pprint(S1)
smat=Matrix(N(S1))

smat=Matrix(N(msq1))
pprint(smat)

print('\nnumerical Svd')
ufloat,sfloat,vfloat = mp.svd(smat)

u=fixfloaterr(ufloat)
#s=fixfloaterr(sfloat)
#v=fixfloaterr(vfloat)

pprint(ufloat)
pprint(u)

pprint(sfloat)
#pprint(s)

pprint(vfloat)
#pprint(v)

print(type(u))
print(type(vfloat))

pprint(Matrix(vfloat))

print()
"""
print('Columns of U \neigenvects of sym*conj(sym)')
pprint((smat*conjugate(smat)).eigenvects())

print('Columns of V\neigenvects of conj(sym)*sym')
pprint((conjugate(smat)*smat).eigenvects())
"""

#unitary=diag(u, eye(2))
#pprint(unitary)
#u1=TensorProduct(eye(2),unitary)

u1=u
v1=Matrix(vfloat)
print('unitary')
pprint(u1)

print('unitary inverse')
pprint(v1)

print('Two mode squeezing matrix')
pprint(msq1)
print('U1 inv * tmsq * u1')
pprint(u1*N(msq1)*v1)





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

