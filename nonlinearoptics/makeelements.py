from sympy import *
from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum import TensorProduct 
from mpmath import mp
from sympy.abc import alpha, beta, xi, zeta
from sympy import I

class Makeelements():
    
    def __init__(self,nspace,nspectral,sq1_mode1,sq1_mode2, sq2_mode1, sq2_mode2, phasemode, phaseangle, bsmode1,bsmode2, bsangle):
   
        self.nspace=nspace
        self.n=nspace*nspectral
        self.nspectral=nspectral
        
        self.a=[None]*self.n
        for i in range(0,self.nspace):
            for j in range(0,self.nspectral):
                self.a[i*(self.nspectral)+j]=symbols('a%d%d' % (i, j))
        
        # phase shift
        omega=[None]*(nspectral)
        for i in range(0,nspectral):
            omega[i]=symbols('omega:%d' % (i))

        omega[0]=pi/2

        ########### Define squeezing symbols
        xi=[None]*(self.n**2)
        for i in range(0,self.n**2):
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

        self.modes=self.makemodes()
      
        #block form for symplectic matrix
        self.al = MatrixSymbol('alpha', self.n,self.n)
        self.be = MatrixSymbol('beta',self.n,self.n)

        #make block form for symplectic matrix
        self.block=BlockMatrix([[self.al,self.be],
                        [conjugate(self.be), conjugate(self.al)]])
 
        #make tmsq on modes 0,1 squeezing xi
        self.sq1=self.makesq(sq1_mode1,sq1_mode2, xi)
        pprint(self.sq1)

        #make tmsq on modes 2,3
        self.sq2=self.makesq(sq2_mode1,sq2_mode2, xi)
        pprint(self.sq2)

        #make phase shifter
        self.ps=self.makeps(phasemode,omega)
        pprint(self.ps)

        #make beamsplitter 
        self.bs=self.makebs(bsmode1,bsmode2,bsangle)
        pprint(self.bs)

        return

    def makemodes(self):
        #Make optical modes
        ablk=MatrixSymbol('a',self.n,1)
        #build block matrix for modes
        bmodes=BlockMatrix([[ablk],[conjugate(ablk)]])
        #make a size n vector for a
        modematrix=Matrix(self.n,1, lambda i,j: self.a[i] )
        #substitution for a modes
        modes=self.matsubs(bmodes,ablk,modematrix,0,0,0,0)
        return modes

    def makeps(self,mode1, phaseangle):
        ##Phase shifter
        m1=mode1*self.nspectral
        phasespace=Matrix(self.n,self.n, lambda i,j: exp(I*phaseangle[i%self.nspectral])
                if ((i==j)and((m1<=i)and(i<m1+self.nspectral))) else 1 if(i==j) else 0)
        # symplectic phase shifer
        mps=Matrix(self.matsubs(self.block,self.al,phasespace,self.be,zeros(self.n),0,0))
        return mps

    def makebs(self,mode1,mode2,theta):
        #Make beamsplitter 
        beamspace= Matrix(self.nspace,self.nspace, lambda i,j: cos(theta) 
                if ((i==j)and((i==mode1)or(i==mode2))) 
                else  sqrt(-1)*sin(theta) if (i+j)==(mode1+mode2)and((i==mode1)or(i==mode2)) 
                else 1 if (i==j)and((i!=mode1)or(i!=mode2)) else 0)
        beamsplitter=TensorProduct(beamspace,eye(self.nspectral))    
        # symplectic beamsplitter
        mbs=Matrix(self.matsubs(self.block,self.al,beamsplitter,self.be,zeros(self.n),0,0))
        return mbs

    def makesinglesq(self):
        #Single mode squeezer
        c1=Matrix(self.n,self.n, lambda i,j: c(xi[i%self.nspectral]) if i==j else 0)
        s1=Matrix(self.n,self.n, lambda i,j: s(xi[i%self.nspectral]) if i==j else 0)
        return c1, s1
        
    def makesq(self,mode1,mode2, sqparam):
        #Two mode squeezer
        #for off-diagonal
        n=self.n
        nspec=self.nspectral
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
        ms01=Matrix(self.matsubs(self.block,self.al,c2,self.be,s2,0,0))
        return ms01

    def matsubs(self,mat,aold,anew,bold,bnew,cold,cnew):
        temp1=mat.subs(aold,anew)
        temp2=temp1.subs(bold,bnew)
        temp3=temp2.subs(cold,cnew)
        return temp3

    def justdoitplease(self,transform, modes, showmodes):
        modetransform=transform*modes
        pprint(relational.Eq(Matrix(modes[0:showmodes]),Matrix(modetransform[0:showmodes])))
        return 

class Matrixops():

    def __init__(self):
        return
    
    def makematrixexact(self,matrix):
        imax,jmax=matrix.shape
        newmat=matrix[:,:] 
        for j in range(0,jmax):
            for i in range(0,imax):
                newmat[i,j]=nsimplify(matrix[i,j], full=True)

        return matrix

    def takagi_for_unitary(self,A):
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

