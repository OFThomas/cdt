
module olis_f90stdlib
!>@file olis_f90stdlib
!> linear alg and random number subroutines
!> oli's standard FORTRAN LIb
!> @author Oliver Thomas
!> August 2018- started docs
!>
!>@brief fortran lib for linear alg
!>
!>@details Allocation routines for temp work arrays for complex_svd,
!> complex_eigenvals and eigenvector solver, outer product of complex vects,
!> matrix formatting print routine and random number and seed subroutines
!>
implicit none

integer, parameter, private :: dp=selected_real_kind(15,300)
real(kind=dp), parameter :: pi=4.0_dp*atan(1.0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! make temp arrays for complex_eigenvects
integer, private :: lda, ldvl, ldvr
integer, parameter, private   :: lwmax = 1000 

! temp scalars
integer, private :: info, lwork

!temp arrays
real(kind=dp), allocatable, dimension(:), private ::  rwork_eigen, rwork_svd
complex(kind=dp), allocatable, dimension(:), private :: work_eigen, work_svd

! make temp arrays for complex_svd

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>@brief allocates eigenvals, u & v arrays for eigenvals & eigenvects
!>
!>@detail allocated temp work arrays also 
!>@author Oliver Thomas
!>August 2018
!>@param matrix input complex matrix 
!>@param eigenvals 1d array for eigenvalues, is overwriten on exit
!>@param u 2d array of left eigenvectors 
!>@param v 3d array of right eigenvectors
subroutine alloc_complex_eigenvects(matrix, eigenvals, u, v)
! complex diag
complex(kind=dp), dimension(:,:), intent(in) :: matrix
complex(kind=dp), dimension(:,:), allocatable, intent(inout) :: u,v
complex(kind=dp), dimension(:), allocatable, intent(inout) :: eigenvals
integer :: n

n=size(matrix,1)

! u, v and eigenvals
allocate(u(n,n))
allocate(v(n,n))
allocate(eigenvals(n))

! temp work arrays
allocate( work_eigen( lwmax ))
allocate(rwork_eigen(2*size(matrix,1)))

print*, 'Allocated temp work arrays for DIAG'
end subroutine alloc_complex_eigenvects

!>@brief allocates sigma (singular vals), u and vt for complexSVD
!>@detail allocates temp work arrays too 
!>@param matrix input complex matrix
!>@param sigma real vector of singular values sorted in descending order
!>@param u unitary matrix
!>@param vt unitary matrix returns V**H NOT v
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine alloc_complex_svd(matrix, sigma, u, vt)
! complex SVD
complex(kind=dp), dimension(:,:), intent(in) :: matrix
complex(kind=dp), dimension(:,:), allocatable, intent(inout) :: u,vt
real(kind=dp), dimension(:), allocatable, intent(inout) :: sigma 
integer :: n,m
m=size(matrix,1)
n=size(matrix,2)

allocate(u(m,n))
allocate(vt( m, n ))
! need to check which is smallest
allocate(sigma(n))

allocate( work_svd( lwmax ))
allocate(rwork_svd(2*size(matrix,1)))

print*, 'Allocated temp work arrays for SVD'
end subroutine alloc_complex_svd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!>@brief generates random seed
!>@param seed is input allocatable 1d array
subroutine randseed(seed)
!!!!! random seed
integer :: values(1:8), seedsize
integer, dimension(:), allocatable :: seed

call date_and_time(values=values)
call random_seed(size=seedsize)

allocate(seed(1:seedsize))

seed(:) = values(8)

call random_seed(put=seed)
end subroutine randseed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11 

!>@brief print formatted matrices
!>@detail can take optional args for labels or write directly to a file
!>@param vect is the input complex matrix
!>@param desc is the optional string to be written above the matrix
!>@param f is the optional file output unit to write to, default is console
subroutine printvectors(vect, desc, f)
  implicit none
  integer :: i, j, m, n
  integer, intent(in), optional :: f 
  character(len=*), intent(in), optional :: desc
  complex(kind=dp), dimension(:,:), intent(in) :: vect
  n=size(vect,1)
  m=size(vect,1)
  if (present(f)) then
          if (present(desc)) then
          write(f,*) desc
          else 
          write(f,*)
          endif 
  do i=1, m
          write(f,9998) (vect (i,j), j=1,n)
  end do
  write(f,*)
  else     
          if (present(desc)) then
          write(*,*) desc
          else 
          write(*,*)
          endif 
  do i=1, m
          write(*,9998) (vect (i,j), j=1,n)
  end do
  write(*,*)
  end if
  
  9998 format( 11(:,1x,'(',F6.2,',',F6.2,')'))
  end subroutine printvectors
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief outerproduct of two complex vectors, returns a complex matrix
  !>@param a is input vector 1, |ket>
  !>@param b is input vector 2, <bra|
  function outerproduct(a,b)
  implicit none
  complex(kind=dp), dimension(2,2) :: outerproduct
  complex(kind=dp), dimension(:), intent(in) :: a, b
  integer :: n, j ,k
  
  ! the return value is the function name
  n=size(a)
  do k=1,n
          do j=1, n
                  outerproduct(j,k)=a(j)*b(k)
          end do
  end do
  end function outerproduct
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !------------------------------ Identity function--------------------------

  !>@brief makes complex identity matrix dim (nxn)
  !>@param n input dimension
  function c_identity(n)
  complex(kind=dp), dimension(n,n) :: c_identity
  integer, intent(in) :: n
  integer :: i
  c_identity=0.0_dp
!make identity matrix nxn
  do i=1,n
    c_identity(i,i)=1
  end do
end function c_identity

!---------------------- Tensor Product function --------------------
!>@brief tensor product for complex matrices aXb
!>@param a complex matrix in
!>@param b complex matrix in
function tprod(a,b)

  complex(kind=dp), dimension (:,:), intent(in) :: a, b
  complex(kind=dp), allocatable, dimension(:,:) :: tprod
  !complex(kind=dp), dimension(:,:) :: tprod
  integer :: ierr, sindex1, sindex2, i,j,k,l, n_a1, n_a2, n_b1, n_b2

  sindex1=size(a, 1)*size(b, 1)
  sindex2=size(a, 2)*size(b, 2)
  
  allocate(tprod(sindex1, sindex2), stat=ierr)
    if (ierr/=0) stop 'Error in allocating tproduct'
  n_a1=size(a,1)
  n_a2=size(a,2)
  n_b1=size(b,1)
  n_b2=size(b,2)

  do j=1, n_a2
    do i=1, n_a1
      do l=1, n_b2
        do k=1, n_b1
          tprod(k+(i-1)*n_b1, l+(j-1)*n_b2) = a(i,j)*b(k,l)
        end do !k
      end do !l
    end do !i
  end do !j

end function tprod

!>@brief computes the trace of a complex matrix
  !>@param a is the complex matrix in
function complextrace(a)
        complex(kind=dp), dimension(:,:) :: a
        complex(kind=dp) :: complextrace 
        integer :: i

        complextrace=0.0_dp
        do i=1, size(a,1)
                complextrace=complextrace+a(i,i)
        end do
end function complextrace
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>@brief computes the complex eigenvalues and eigenvectors
!>@detail overwrites matrix in, input eigenvalue array and eigenvector arrays
!> uses the zgeev subroutine from lapack
!>@param a input allocatable complex  matrix to be diagonalised
!>@param w output allocatable complex 1d array containing eigenvals
!>@param vl output allocatable complex 2d array containing left eigenvectors
!>@param vr output allocatable complex 2d array containing right eigenvectors
!>@note need to check this is optimised
subroutine complex_eigenvects(a, w, vl, vr)
! matrix a in, eigenvals, eigenvects out
implicit none
! matrix in is a
!!!!!!!!!!!!!!!!!!!!!!! matrix a is overwritten WATCH OUT!
complex(kind=dp), dimension(:,:), allocatable :: a

! left & right vectors
complex(kind=dp), dimension(:,:), allocatable :: vl
complex(kind=dp), dimension(:,:), allocatable :: vr
! eigen values are w
complex(kind=dp), allocatable, dimension(:) :: w 

integer :: n, m

! use size of input matrix
n=size(a,1)
m=size(a,2)

lda=n
ldvl=n
ldvr=m

!     qUERY THE OPTIMAL WORKSPACE.
lwork = -1
call zgeev( 'V', 'v', n, a, lda, w, vl, ldvl, vr, ldvr, work_eigen, lwork, rwork_eigen, info )
lwork = min( lwmax, int( work_eigen( 1 ) ) )

!     sOLVE EIGENPROBLEM.
call zgeev( 'v', 'v', n, a, lda, w, vl, ldvl, vr, ldvr, work_eigen, lwork, rwork_eigen, info )

!     cHECK FOR CONVERGENCE.
if( info.gt.0 ) then
write(*,*)'tHE ALGORITHM FAILED TO COMPUTE EIGENVALUES.'
stop
end if
!if (n==2) then
!    vr=c_inv2(vl)
!end if
end subroutine complex_eigenvects

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>@brief computes the complex eigenvalues and eigenvectors
!>@detail overwrites matrix in, input eigenvalue array and eigenvector arrays
!> uses the zgeev subroutine from lapack
!>@param a input allocatable complex  matrix to be SVD'd
!>@param sigma output allocatable complex 1d array containing ordered singular values
!>@param u output allocatable complex 2d array containing u
!>@param vt output allocatable complex 2d array containing v**H
!>@note need to check this is optimised
subroutine complex_svd(a, sigma, u, vt)
! matrix a in, eigenvals, eigenvects out
! A = U * sigma * V **H
implicit none 

integer :: n, m, i
integer :: lda, ldu, ldvt
integer, parameter   :: lwmax = 1000 

! matrix in is a
complex(kind=dp), dimension(:,:), allocatable, intent(inout) :: a
! singular values sigma
real(kind=dp), dimension(:), allocatable :: sigma
! left & right vectors
complex(kind=dp), dimension(:,:), allocatable :: u, vt
! eigen values are w

! use size of input matrix
! might have got n & m the wrong way round
n=size(a,1)
m=size(a,2)
ldu=size(a,1)
ldvt=size(a,1)
lda=size(a,1)

! eigen vectors, eigen vals & temp arrays

!no left and right col vectors
! rows m =size(a,1) ! cols n= size(a,2) ! a is matrix ! lda =size(a,1)
! s vector svd ! u matrix ! ldu = m  ! vt matrix hermitian conjg
! ldvt = n  ! work ! lwork! rwork
call printvectors(a, 'dens -dens_est')

! quiery the workspace size
lwork = -1
call zgesvd('S','S', m, n, a, lda, sigma, u, ldu, vt, ldvt, work_svd, lwork, rwork_svd, info )

! do svd
lwork = min( lwmax, int( work_svd( 1 ) ) )
call zgesvd('S','S', m, n, a, lda, sigma, u, ldu, vt, ldvt,  work_svd, lwork, rwork_svd, info )

!     cHECK FOR CONVERGENCE.
if( info.gt.0 ) then
write(*,*)'tHE ALGORITHM FAILED TO COMPUTE EIGENVALUES.'
stop
end if

! write singular vals back to matrix a
a=0.0_dp
do i =1, n
    a(i,i)=sigma(i)
end do
end subroutine complex_svd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>@brief inverse for a complex 2x2 matrix
!>@param m_in is input complex 2x2 matrix
function c_inv2(m_in)
    implicit none
    complex(kind=dp), dimension(2,2), intent(in) :: m_in 
    complex(kind=dp), dimension(2,2) :: c_inv2
    complex(kind=dp) :: det
    c_inv2(1,1)=m_in(2,2)
    c_inv2(1,2)=-m_in(1,2)
    c_inv2(2,1)=-m_in(2,1)
    c_inv2(2,2)=m_in(1,1)
    det=m_in(1,1)*m_in(2,2) - m_in(1,2)*m_in(2,1)
    c_inv2=c_inv2/det
end function c_inv2


!>@brief computed Frobenieus matrix norm of complex matrix using lapack zlange
!>@param c input complex matrix 
function matrixnorm(c)
    complex(kind=dp), dimension(:,:) :: c
    real(kind=dp) :: matrixnorm, zlange
    ! temp
    real(kind=dp), dimension(:), allocatable :: work
    integer :: m,n,lda, lwmax=1000 
 
    m=size(c,1)
    n=size(c,2)
    lda=m
    allocate(work(lwmax))
    matrixnorm= zlange('F', m,n,c,lda,work )
end function matrixnorm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
end module olis_f90stdlib
