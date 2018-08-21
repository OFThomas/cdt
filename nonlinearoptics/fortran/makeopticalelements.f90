!>@brief module for building symplectic matrices for optical elements
module makeopticalelements
!>@file makeopticalelements 
!> matrix builers and memory allocation for optical symplectic matrices
!>
!> @author Oliver Thomas
!> August 2018- started docs
!>
use olis_f90stdlib 
implicit none

integer, parameter, private :: dp=selected_real_kind(15,300)
real(kind=dp), public :: ident

!>@param temp work arrays for beamsplitter 
complex(kind=dp), dimension(:,:), allocatable, private :: ident_spec, spatial_work, n_work

contains 
    !>@brief makes beamsplitter symplectic matrix
    !>@detail takes in an allocated matrix for the beamsplitter matrix
    !> to be written to
    !> uses the private ident_spec, spatial_work, n_work arrays 
    !>@param nspace is number of total spatial modes
    !>@param nspec is number of total spectral modes
    !>@param m_bs allocated n*n matrix for beamsplitter
    !>@param m1 is spatial mode 1 for beam splitter
    !>@param m2 is spatial mode 2 for beam splitter
    subroutine make_bs(nspace,nspec, symp_mat, m1, m2, theta)
        integer :: nspace, nspec, n, i, j,m1,m2
        real(kind=dp) :: theta
        ! temp work matrices
        !complex(kind=dp), dimension(:,:), allocatable :: ident_spec, spatial_work, n_work

        complex(kind=dp), dimension(:,:), allocatable :: symp_mat
        
        n=nspace*nspec
        !allocate(m_bs(nspace,nspace))
        spatial_work=0.0_dp
        do j=1, nspace
            ! set all diag to ident
            spatial_work(j,j)=1.0_dp
            do i=1, nspace
                ! if in desired spatial mode do 
                if ((i == m1) .or. (i == m2)) then
                    ! on diag is cos(theta)
                    if (i==j) then
                        spatial_work(i,j)=cos(theta)
                    ! off diag is i sin(theta)
                    elseif ((i + j) == (m1 + m2)) then 
                        spatial_work(i,j)=complex(0.0_dp,sin(theta))
                    end if
                end if
            end do
        end do
    ident_spec=c_identity(nspec)
    ! do tensor product over all spectral modes
    n_work=tprod(spatial_work,ident_spec)
    ! returns alpha
    ! now to get symplectic matrix
    symp_mat=0.0_dp
    symp_mat(1:n,1:n)=n_work
    symp_mat(n+1:2*n,n+1:2*n)=conjg(n_work)
end subroutine make_bs

!>@brief make symplectic squeezing matrix from exponetiated JSA
!>@TODO a lot is broken...
!>@note only works if modes are consectutive
subroutine make_sq(nspace,nspec,symp_mat,m1,m2, alpha, beta)

    integer :: nspace, nspec, m1, m2
    complex(kind=dp), dimension(:,:), allocatable :: symp_mat, a, b
    complex(kind=dp), dimension(:,:), intent(inout) :: alpha, beta
    
    !counters
    integer :: i,j,k,l, n
    integer :: m1s, m2s

    n=nspace*nspec
    m1s=m1*nspec
    m2s=m2*nspec
    a=alpha
    b=beta

    !if ((m1s <= i < m1s + nspec) .or. (m2s <= i < m2s + nspec)) then
    !>@note alpha & beta are 2 spatial modes and all spectral modes
    !> dim 2*nspace*nspec
   
    ! make passives identity 
    symp_mat(:,:)=c_identity(n*2)

    !!then overwrite values from alpha
    !! do alpha first on modes m1, m2 
    !! for m1 do for all spectral modes copy alpha
    
    !do j=m1s,m1s+nspec
        !do i=m1s,m1s+nspec-1
            
            !symp_mat(i,j)=alpha(i-m1s+1,j-m1s+1)
            !print*, 'i', i, 'j', j
            !print*, 'i-m1s+1', i-m1s+1, 'j-m1s+1', j-m1s+1
        !end do
    !end do
    
    !do j=m2s,m2s+nspec
    !    do i=m2s, m2s+nspec-1
    !        symp_mat(i,j)=alpha(i-m2s+1,j-m2s)
    !        print*, 'i2', i, 'j2', j
    !        print*, 'i-m1s+1', i-m1s+1, 'j-m1s', j-m1s
    !    end do
    !end do
    
    !>@brief loop for alpha 
    l=1
    do j=m1s, m2s
        k=1
        do i=m1s, m2s
            symp_mat(i,j)=alpha(k,l)
            k=k+1
        end do
        l=l+1
    end do

    !symp_mat(m1s:m1s+nspec, m1s:m1s+nspec)=alpha(1:nspec,1:nspec)
    !1:nspec)
    ! for m2 do the same
    !symp_mat(m2s:m2s+nspec, m2s:m2s+nspec)=alpha(nspec+1:2*nspec, nspec+1:2*nspec)
    
    
    ! alpha conjg is then 
    symp_mat(n+1:2*n, n+1:2*n)=conjg(symp_mat(1:n, 1:n))

    ! beta is much easier
    !same as alpha but second index is off-set by n
   

    ! for mode 1
    !>@TODO check this is legal...
    !> full diag sq
    !> symp_mat(m1s:m1s+nspec, m1s+n:m1s+nspec+n)=beta(1:nspec, 1+nspec:2*nspec)
   
    ! mode 1
    !symp_mat(m1s:m1s+nspec, m2s+n:m2s+nspec+n)=beta(1:nspec, 1+nspec:2*nspec)
    
    ! for mode 2
    !symp_mat(m2s:m2s+nspec, m1s+n:m1s+nspec+n)=beta(nspec+1:2*nspec, 1:nspec)
    
    !>@TODO
    !> probably not legal
    !>symp_mat(m2s:m2s+nspec, m2s+n:m2s+nspec+n)=beta(nspec+1:2*nspec, 1:nspec)


    !>@brief loop for beta, offset to col+n
    l=1
    do j=m1s+n,m2s+n
        k=1
        do i=m1s,m2s
            symp_mat(i,j)=beta(k,l)
            k=k+1    
        end do
    l=l+1
    end do
    ! beta conjg is then
    symp_mat(n+1:2*n, 1:n)=conjg(symp_mat(1:n, n+1:2*n))


end subroutine make_sq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>@brief make sqq matrix from jsa function
function make_squeezer(nspace, nspec, mode1, mode2, jsa) 

integer, intent(in) :: mode1, mode2, nspace, nspec
complex(kind=dp), dimension(:,:), allocatable, intent(in) :: jsa
complex(kind=dp), dimension(:,:), allocatable :: make_squeezer


complex(kind=dp), dimension(:,:), allocatable :: h_sq, alpha_temp, beta_temp
integer :: f_size, sizeofexp, alpha_size, beta_size
f_size=size(jsa,1)
sizeofexp=4*f_size
allocate(h_sq(sizeofexp,sizeofexp))

!!!!!!!!!
!>@note to make off diagonal for fmatrix
!>m_sq=0.0_dp
!>! top right
!>m_sq(1:1*f_size, 3*f_size+1:4*f_size)=1
!>! mid right
!>m_sq(1*f_size+1:2*f_size, 2*f_size+1:3*f_size)=2 
!>! mid left
!>m_sq(2*f_size+1:3*f_size, 1*f_size+1:2*f_size)=3
!>! bot left
!>m_sq(3*f_size+1:4*f_size, 1:1*f_size)=4
!call printvectors(m_sq)

!>@note !h= 0.0       F_JSA
!>           F_JSA*T   0.0
!>
!>f_jsa = f_mat
!>
!> M_sq = exp(i ( 0 H )
!>              (-H* 0)
!>
!> M_sq = exp(i  (0              0           0       F_JSA) 
!>               (0              0         F_JSA**T      0)
!>               (0         -conjg(F_JSA)    0           0)
!>               (-F_JSA**H      0           0           0)

h_sq=0.0_dp
! top right
h_sq(1 : 1*f_size,            (3*f_size)+1 : 4*f_size)=jsa(:,:)
! mid right
h_sq((1*f_size)+1 : 2*f_size, (2*f_size)+1 : 3*f_size)=transpose(jsa(:,:)) 
! mid left
h_sq((2*f_size)+1 : 3*f_size, (1*f_size)+1 : 2*f_size)=-conjg(jsa(:,:))
! bot left
h_sq((3*f_size)+1 : 4*f_size,  1 : 1*f_size)=-conjg(transpose(jsa(:,:)))
!call printvectors(m_sq)

! do exponentiation
h_sq=expmatrix(imaginary*h_sq,50)
!call printvectors(m_sq, 'after exp')

alpha_size=2*f_size
allocate(alpha_temp(alpha_size,alpha_size))
allocate(beta_temp(alpha_size,alpha_size))

!>@note alpha beta are top left and top right of M
!> M = (A  B )
!>     (B* A*)
!>@param alpha_size is 2*f_size as all spectral modes for 2 spatial
alpha_temp(:,:)=h_sq(1:alpha_size,1:alpha_size)
beta_temp(:,:)=h_sq(1:alpha_size, alpha_size+1:2*alpha_size)

!call printvectors(alpha_temp, 'alpa')
!call printvectors(beta_temp, 'beta')

!@brief inset squeezing hamiltonian in correct spatial modes
!>@note allocate for sq on modes 1&2
allocate(make_squeezer(2*nspace*f_size,2*nspace*f_size))
make_squeezer=0.0_dp
! 0.25 of matrix, need other 0.25 to be ident
! half of the alpha block

! check how many modes there are,
! if more modes than alpha make the rest of sq ident

! do for modes 3&4
! 0 to 2*n

!allocate(mat_sq2(2*nspace*f_size, 2*nspace*f_size))
!mat_sq2=0.0_dp


call make_sq(nspace, f_size, make_squeezer,mode1,mode2,alpha_temp,beta_temp)
!call printvectors(mat_sq1, 'sq on modes 1&2')

!call printvectors(alpha_temp, 'alpa')
!call printvectors(beta_temp, 'beta')

!call make_sq(nspace, f_size, mat_sq2, 3,4,alpha_temp,beta_temp)
!call printvectors(mat_sq2, 'sq on modes 3&4')

!call printvectors(alpha_temp, 'alpa')
!call printvectors(beta_temp, 'beta')

end function make_squeezer


!>@brief samples the given jsa for frequency ranges w1, w2
!>@param f_mat allocatable Jsa matrix values out
!>@param w_start, 
!>@param w_end
!>@param w_incr
!>@param sigma is jsa parameter
!>@param outfile is unit number of file to write to 
function gen_jsa(f,w1_start, w1_end, w1_incr, w2_start, w2_end, w2_incr, sigma1, sigma2, outfile, w1offset, w2offset)
real(kind=dp), dimension (:,:), allocatable :: gen_jsa
real(kind=dp), intent(in) :: w1_start, w1_end, w1_incr
real(kind=dp), intent(in) :: w2_start, w2_end, w2_incr
real(kind=dp), intent(in) :: sigma1, sigma2 
real(kind=dp) :: w1, w2
real(kind=dp) :: w1off, w2off
integer :: w1_steps, w2_steps, outfile
integer :: i, j

complex(kind=dp) :: f


real(kind=dp), optional :: w1offset, w2offset

w1off=0.0_dp; w2off=0.0_dp

if(present(w1offset)) w1off=w1offset
if(present(w2offset)) w2off=w2offset

w1_steps=ceiling((w1_end-w1_start)/w1_incr)
w2_steps=ceiling((w2_end-w2_start)/w2_incr)
print*, 'dim of JSA is', w1_steps, w2_steps

if (w1_steps >= w2_steps) then 
    allocate(gen_jsa(w1_steps,w1_steps))
else if (w1_steps < w2_steps) then
    allocate(gen_jsa(w2_steps,w2_steps))
else 
    stop 'what the hell is wrong with the spectral range'
end if

!do l=1,2
w2=w2_start
do j=1,w2_steps
    w1=w1_start
    do i= 1,w1_steps 
        write(outfile,*) w1, w2, real(f(w1,w2,sigma1, sigma2, w1off,w2off),kind=dp)
        gen_jsa(i,j)=real(f(w1,w2,sigma1, sigma2, w1off, w2off),kind=dp)
        w1=w1+w1_incr
    end do
    w2=w2+w2_incr
end do
end function gen_jsa 

!>@brief JSA function taking two freq
!>@param w1 input signal freq
!>@param w2 input idler freq
!>@param sig input variance
    function f_gauss(w1,w2, sigma1, sigma2, w1off, w2off)
    complex(kind=dp) :: f_gauss
    real(kind=dp), intent(in) :: w1,w2, sigma1, sigma2, w1off, w2off
    
    f_gauss=(1.0_dp/(sigma1*sigma2)) * exp(-0.5_dp*((w1-w1off)/sigma1)**2)*exp(-0.5_dp*((w2-w2off)/sigma2)**2)

    end function f_gauss

    function f_sine(w1,w2,sigma1,sigma2,w1off,w2off)
    complex(kind=dp) :: f_sine
    real(kind=dp), intent(in) :: w1,w2,sigma1,sigma2, w1off,w2off

    f_sine=sinc(2.0_dp*(w1-w2))*exp(-0.5*(w1**2+(w2*2.0_dp)**2))
    end function f_sine

!>@brief calculates g4 using matrix elements sum
!>@TODO
!>@param ft input is the full symplectic transform
!>@param nspec input spectral DOF
function g4(ft, nspec)
real(kind=dp) :: g4
complex(kind=dp), dimension(:,:), allocatable, intent(in) :: ft
integer, intent(in) :: nspec
complex(kind=dp), dimension(6) :: term
complex(kind=dp) :: gam21, gam43, gam32, gam41, gam31, gam42
real(kind=dp) :: bbdag11, bbdag22, bbdag33, bbdag44
real(kind=dp) :: gam, bdiag
g4=0.0_dp

!print*, 'did abt'
gam21 = abt(2,1,ft, nspec)
gam43 = abt(4,3,ft, nspec)
gam32 = abt(3,2,ft, nspec)
gam41 = abt(4,1,ft, nspec)
gam31 = abt(3,1,ft, nspec)
gam42 = abt(4,2,ft, nspec)

gam=amp(gam21*gam43 + gam32*gam41 + gam31*gam42)

bbdag11 = bbd(1,1,ft, nspec)
bbdag22 = bbd(2,2,ft, nspec)
bbdag33 = bbd(3,3,ft, nspec)
bbdag44 = bbd(4,4,ft, nspec)

!print*, 'gam21 abt'
bdiag= bbdag11*bbdag22*bbdag33*bbdag44

term(1) = amp(gam21) * bbdag44 * bbdag33
term(2) = amp(gam32) * bbdag11 * bbdag44
term(3) = amp(gam31) * bbdag22 * bbdag44
term(4) = amp(gam43) * bbdag11 * bbdag22
term(5) = amp(gam42) * bbdag11 * bbdag33
term(6) = amp(gam41) * bbdag22 * bbdag33

g4 = abs(gam + bdiag + term(1) + term(2) + term(3) + term(4) + term(5) +term(6))
print*, 'gam', gam, 'bdiag', bdiag
end function g4

!>@brief returns the absolute value squared |a|**2
!>@param a input complex number to be |a|**2
function amp(a)
real(kind=dp) :: amp
complex(kind=dp) ::a
amp=abs(a)**2
end function amp

!>@brief calculates matrix elements Alpha-Beta**T
!>@detail for M = (A  B )
!>                (B* A*)
!> computes AB**T and returns the i,j-th element
!>@param i input index 1
!>@param j input index 2
!>@param ft input symplectic transform matrix for the optical circuit
!>@param nspec input number of spectral DOF
function abt(i,j,ft, nspec)
complex(kind=dp) :: abt
integer :: i,j, nspec
integer :: n, a, b, k
complex(kind=dp), dimension(:,:), allocatable, intent(in) :: ft
complex(kind=dp), dimension(:,:), allocatable :: temp
!total dim 2*n
n=nint(0.5_dp*size(ft,1))

if ((i < n).and.(j<n)) then 
    a=1
    b=n
else if ((i>=n).and.(j>=n)) then
    a=n+1
    b=2*n
    i=i-n
    j=j-n
end if

temp=matmul( ft(a:b, 1:n), transpose(ft(a:b, n+1:2*n)) ) 
abt=0.0_dp
!>@todo check this
do k=1,nspec
    abt=abt+temp(((i-1)*nspec)+k, ((j-1)*nspec)+k)
end do
!abt=abt/real(nspec,kind=dp)
end function abt


!>@brief calculates the matrix elements Beta*Beta**H
!>@detail for M = (A  B )
!>                (B* A*)
!> computes B*B**H (Hermitian conjg) and returns the i,j-th element
!>@param i input index 1
!>@param j input index 2
!>@param ft input symplectic transform matrix for the optical circuit
!>@param nspec input number of spectral DOF
function bbd(i,j,ft,nspec)
complex(kind=dp) :: bbd
integer, intent(in) :: i,j, nspec
integer :: n, k
complex(kind=dp), dimension(:,:), allocatable, intent(in) :: ft
complex(kind=dp), dimension(:,:), allocatable :: temp

n=nint(0.5_dp*size(ft,1))
!print*, 'n', n, 'nspec', nspec
!print*, 'i', i, 'j', j
temp=matmul(ft(1:n, n+1: 2*n), conjg(transpose(ft(1:n, n+1:2*n))))
bbd=0.0_dp
!>@todo check this
do k=1, nspec
    bbd=bbd+temp(((i-1)*nspec)+k, ((j-1)*nspec)+k)
    !print*, 'bbd sum', bbd
    !print*, ((i-1)*nspec)+k, ((j-1)*nspec)+k
end do
!bbd=bbd/real(nspec,kind=dp)
end function bbd






!>@brief allocates temp arrays for matrices
!>@param nspace input
!>@param nspec input
!>@detail allocates memory for ident_spec a spectral size matrix for tensor producting.
!>
!>allocates mem for spatial_work, array size of spatial modes
!>
!> allocates mem for n_work, work array size of alpha or beta in sympectic matrix
subroutine alloc_temparrays(nspace,nspec)
integer, intent(in) :: nspace, nspec

allocate(ident_spec(nspec,nspec))
allocate(spatial_work(nspace,nspace))
allocate(n_work(nspace*nspec,nspace*nspec))

end subroutine alloc_temparrays

subroutine dealloc_temparrays

deallocate(ident_spec)
deallocate(spatial_work)
deallocate(n_work)
end subroutine dealloc_temparrays




end module makeopticalelements
