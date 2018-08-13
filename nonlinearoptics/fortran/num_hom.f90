!>@brief program to compute matrix of a JSA
program num_hom
use olis_f90stdlib
use makeopticalelements 
implicit none

!>@param dp 
integer, parameter :: dp=selected_real_kind(15,300)

!>@param i counter
!>@param j counter
!>@param k counter
integer :: i,j,k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! params for plotting f(w1,w2)
integer :: w1_steps, w2_steps

!>@param w1 signal freq
!>@param w2 idler freq
!>@param sigma is gaussian width
real(kind=dp) :: w1, w2, sigma

! loop params
!>@param w1_start the start frequency range for signaler
!>@param w2_start the start frequency range for idler
!>@param w1_incr spectral dof mesh
!>@param w2_incr spectral dof mesh
real(kind=dp) :: w1_start, w1_end, w2_start, w2_end
real(kind=dp) :: w1_incr, w2_incr
!>@param f_mat matrix for values of function, f
!>@param mat_bs symplectic matrix for the beamsplitter
!>@param mat_sq1 symplectic matrix for the squeezer on mode 1&2
!>@param mat_sq2 syplectic matrix for squeezer on modes 3&4
complex(kind=dp), dimension(:,:), allocatable :: f_mat, mat_bs, mat_sq1, mat_sq2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>@param nspace number of spatial modes
!>@param nspec number of spectral modes per spatial mode
!>@param n=nspace*nspec = total number of all modes
integer :: nspace, nspec, n


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! for JSA calcs
!>@param h 
!>@param m_sq 
!>@param m_sq1
!<@param m_sq2 h 
!>@param m_sq 
!>@param m_sq1
!<@param m_sq2
complex(kind=dp), dimension(:,:), allocatable :: h, m_sq, m_sq1, m_sq2
integer :: sizeofexp, f_size, alpha_size

!>@param temp vaaibles for building exp on correct spatial modes
complex(kind=dp), dimension(:,:), allocatable :: alpha_temp, beta_temp

!>@param transform allocatable array for storing total symplectic matrix transform in
!> temp for transform
!>@param transform_no_bs array for storing symplectic transform without beam splitter for normalising g4  
complex(kind=dp), dimension(:,:), allocatable :: transform, transform_no_bs

!>@param theta beamsplitter angle
!>@param g4norm is the normalised g4 matrix
real(kind=dp) :: theta, g4norm, g4_nobs, g4un_norm



!!!!
complex(kind=dp), dimension(:,:), allocatable :: f_mat1, f_mat2
!>@note files to write to  
open(unit=14,file='fplotw1w2.dat', status='replace')
open(unit=15,file='fplotw3w4.dat', status='replace')
open(unit=16, file='g4f90data.dat', status='replace')

w1_start=-5.0_dp
w2_start=-5.0_dp

w1_end=-w1_start
w2_end=-w2_start

!0.05 
w1_incr=0.15_dp
w2_incr=0.15_dp

!w1_steps=ceiling((w1_end-w1_start)/w1_incr)
!w2_steps=ceiling((w2_end-w2_start)/w2_incr)

sigma=1.0_dp

!allocate(f_mat(w1_steps,w2_steps))

!do l=1,2

!w2=w2_start
!do j=1,w2_steps
!    w1=w1_start
!    do i= 1,w1_steps 
!        write(15,*) w1, w2, real(f(w1,w2,sigma),kind=dp)
!        f_mat(i,j)=real(f(w1,w2,sigma),kind=dp)
!        w1=w1+w1_incr
!    end do
!    w2=w2+w2_incr
!end do


f_mat1= gen_jsa(w1_start, w1_end, w1_incr, w2_start, w2_end, w2_incr, sigma,14)
f_mat2=gen_jsa(w1_start,w1_end,w1_incr,w2_start,w2_end,w2_incr,2.0_dp,15)

! normalise?
f_mat1=f_mat1/sum(f_mat1)
f_mat2=f_mat2/sum(f_mat2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 




nspace=4
nspec=size(f_mat1,1)
print*, 'nspectral dof', nspec
n=nspace*nspec
allocate(mat_bs(2*n,2*n))

call alloc_temparrays(nspace,nspec)

!call printvectors( matmul(mat_bs,conjg(transpose(mat_bs))), 'print bs pi/4')
!do i=1, size(f_mat,1)
!    write(15,*) j,i, (real(f_mat(i,j)), j=1, size(f_mat,1))
!end do

f_size=size(f_mat1,1)
sizeofexp=4*f_size
print*, 'dim of exp', sizeofexp
print*, 'transform_no_bs'


!
mat_sq1=make_squeezer(1,2,f_mat1)

!
mat_sq2=make_squeezer(3,4,f_mat2)


transform_no_bs=matmul(mat_sq1,mat_sq2)
print*, 'nspec',  nspec
g4_nobs=g4(transform_no_bs,nspec)

theta=0.0_dp
do i=1, 200
    theta=theta+0.01_dp*pi

    ! modes 2 & 3 
    call make_bs(nspace,nspec,mat_bs,2,3,theta)
    !call printvectors(mat_bs, 'beamsplitter')


!transform_no_bs=matmul(mat_sq1,mat_sq2)
!g4_nobs=g4(transform_no_bs,nspec)

    !call printvectors(transform, 'sq1 * sq2')

    !print*, 'now do beam splitter on modes 2&3'

    transform=matmul(mat_bs,transform_no_bs)
    !call printvectors(transform(1:n,1:n), 'sq1*sq2*BS') 
    g4un_norm= g4(transform, nspec)
    g4norm=g4un_norm/g4_nobs
    print*, 'theta', theta, 'g4 no BS', g4_nobs, 'g4 un',g4un_norm, 'g4norm', g4norm
    write(16,*) theta, g4_nobs,g4un_norm

end do

deallocate(mat_bs)
!deallocate(alpha_temp)
!deallocate(beta_temp)
deallocate(mat_sq1)
deallocate(mat_sq2)
call dealloc_temparrays()

!end do
close(14)
close(15)
close(16)
contains 

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>@brief make sqq matrix from jsa function
function make_squeezer(mode1, mode2, jsa) 

integer, intent(in) :: mode1, mode2
complex(kind=dp), dimension(:,:), allocatable, intent(in) :: jsa
complex(kind=dp), dimension(:,:), allocatable :: make_squeezer

complex(kind=dp), dimension(:,:), allocatable :: h_sq, alpha_temp, beta_temp
integer :: f_size, sizeofexp
f_size=size(f_mat1,1)
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
h_sq(1:1*f_size, 3*f_size+1:4*f_size)=jsa(:,:)
! mid right
h_sq(1*f_size+1:2*f_size, 2*f_size+1:3*f_size)=transpose(jsa(:,:)) 
! mid left
h_sq(2*f_size+1:3*f_size, 1*f_size+1:2*f_size)=-conjg(jsa(:,:))
! bot left
h_sq(3*f_size+1:4*f_size, 1:1*f_size)=-conjg(transpose(jsa(:,:)))
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11


!>@brief samples the given jsa for frequency ranges w1, w2
!>@param f_mat allocatable Jsa matrix values out
!>@param w_start, 
!>@param w_end
!>@param w_incr
!>@param sigma is jsa parameter
!>@param outfile is unit number of file to write to 
function gen_jsa(w1_start, w1_end, w1_incr, w2_start, w2_end, w2_incr, sigma, outfile)
real(kind=dp), dimension (:,:), allocatable :: gen_jsa
real(kind=dp), intent(in) :: w1_start, w1_end, w1_incr
real(kind=dp), intent(in) :: w2_start, w2_end, w2_incr
real(kind=dp), intent(in) :: sigma 
real(kind=dp) :: w1, w2
integer :: w1_steps, w2_steps, outfile
integer :: i, j

w1_steps=ceiling((w1_end-w1_start)/w1_incr)
w2_steps=ceiling((w2_end-w2_start)/w2_incr)
allocate(gen_jsa(w1_steps,w2_steps))
!do l=1,2
w2=w2_start
do j=1,w2_steps
    w1=w1_start
    do i= 1,w1_steps 
        write(outfile,*) w1, w2, real(f(w1,w2,sigma),kind=dp)
        gen_jsa(i,j)=real(f(w1,w2,sigma),kind=dp)
        w1=w1+w1_incr
    end do
    w2=w2+w2_incr
end do
end function gen_jsa 

!>@brief JSA function taking two freq
!>@param w1 input signal freq
!>@param w2 input idler freq
!>@param sig input variance
    function f(w1,w2, sig)
    complex(kind=dp) :: f
    real(kind=dp), intent(in) :: w1,w2, sig
    
    f=(3.0_dp/(sig*sqrt(2.0_dp*pi)))**2 * exp(-0.5_dp*(w1/sig)**2)*exp(-0.5_dp*(w2/sig)**2)

    end function f
end program num_hom
