!>@brief program to compute matrix of a JSA
program num_hom
use olis_f90stdlib
use makeopticalelements 
use schmidt_decomp
implicit none

!>@param dp 
integer, parameter :: dp=selected_real_kind(15,300)

!>@param i counter
!>@param j counter
!>@param k counter
integer :: i,j,k, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! params for plotting f(w1,w2)

!>@param sigma1 is gaussian width for signal
!>@param sigma2 is gaussian width for idler
real(kind=dp) :: sigma1, sigma2

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
complex(kind=dp), dimension(:,:), allocatable :: mat_bs, mat_sq1, mat_sq2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>@param nspace number of spatial modes
!>@param nspec number of spectral modes per spatial mode
!>@param n=nspace*nspec = total number of all modes
integer :: nspace, n
! 2 spectral photons signal #& idler
integer, dimension(2) :: nspec


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! for JSA calcs
!>@param m_sq list of arrays of squeezing matrices
complex(kind=dp), dimension(:,:,:), allocatable :: m_sq

!>@param sizeofexp is size of the whole symplectic space
!>@param f_size is number of dicrete freq vals
!>@param alpha_size is the size of the passive block in symplectic matrix
integer :: sizeofexp, f_size, alpha_size

!>@param temp variables for building exp on correct spatial modes

!>@param transform allocatable array for storing total symplectic matrix transform in
!> temp for transform
!>@param transform_no_bs array for storing symplectic transform without beam splitter for normalising g4  
complex(kind=dp), dimension(:,:), allocatable :: transform, transform_no_bs

!>@param theta beamsplitter angle
!>@param g4norm is the normalised g4 matrix
real(kind=dp) :: thetaincr, g4_nobs, g4un_norm, jsaamp
real(kind=dp), dimension(:), allocatable :: g4norm, theta

integer :: theta_samples


!!!!
complex(kind=dp), dimension(:,:), allocatable :: f_mat1, f_mat2
complex(kind=dp), dimension(:,:,:), allocatable :: f_mat
!>@param uf1 is u matrix from f_mat1 svd
!>@param vtf1 is vt matrix from f_mat1 svd
complex(kind=dp), dimension(:,:,:), allocatable :: uf, vtf, uf2, vtf2
!>@param svf1 singular values for f_mat1
real(kind=dp), dimension(:,:), allocatable :: svf, svf2

!>@param signalfreq array of schmidt decomp values (u matrix)
real(kind=dp) :: signalfreq
!>@param idlerfreq array of schmidt decomp values (vt matrix)
real(kind=dp) :: idlerfreq


integer :: f_dim1, f_dim2 ,num_sq
integer, dimension(:,:), allocatable :: units



!>@note files to write to  
open(unit=14,file='fplotw1w2.dat', status='replace')
open(unit=15,file='fplotw3w4.dat', status='replace')

!open(unit=16, file='g4f90data.dat', status='replace')
!open(unit=17, file='g4splot.dat', status='replace')
!>
!> jsa 1 
open(unit=20, file='signalfreq1.dat', status='replace')
open(unit=21, file='idlerfreq1.dat', status='replace')
!> jsa 2
open(unit=22, file='signalfreq2.dat', status='replace')
open(unit=23, file='idlerfreq2.dat', status='replace')

num_sq=2

w1_start=-5.0_dp
w2_start=-5.0_dp

w1_end=-w1_start
w2_end=-w2_start

!0.05 
w1_incr=0.10_dp
w2_incr=0.10_dp

sigma1=1.0_dp
sigma2=0.5_dp*sigma1


f_dim1=nint((w1_end-w1_start)/w1_incr)+1
f_dim2=nint((w2_end-w2_start)/w2_incr)+1
allocate(f_mat(f_dim1,f_dim2,num_sq))

f_mat(:,:,1)= gen_jsa(f_sine, w1_start, f_dim1, w1_incr, w2_start, f_dim2, w2_incr, &
    sigma1, sigma2,14, w1offset=0.0_dp, w2offset=0.0_dp)

f_mat(:,:,2)=gen_jsa(f_sine, w1_start, f_dim1, w1_incr, w2_start, f_dim2, w2_incr, &
    sigma1, sigma2,15, w1offset=1.0_dp, w2offset=1.0_dp)

! normalise?
do i =1, num_sq
f_mat(:,:,i)=f_mat(:,:,i)/(sum(f_mat(:,:,i))*w1_incr*w2_incr)
write(*,*) 'sum f1',i,  sum(f_mat(:,:,i))*w1_incr*w2_incr
print*, 's1 f', size(f_mat,1), 'f_dim', f_dim1, f_dim2
print*, 'total f size', size(f_mat)
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

nspace=4
nspec(1)=f_dim1
nspec(2)=f_dim2
print*, 'nspectral dof', nspec
n=nspace*maxval(nspec)


! jsa to be decomp, singular vals,
! u, vt
!>@brief allocates the singular values, u and vt matrices for svd

call alloc_temparrays(nspace,nspec(1))

allocate(units(2,num_sq))
units(:,1)=(/20,21/)
units(:,2)=(/22,23/)

!>@note call schmidt_modes(jsa_func, sv, u, vt, units )

!>@note this is wrong...
allocate(svf(nspec(1),num_sq))
allocate(uf(nspec(2),nspec(2),num_sq))
allocate(vtf(nspec(2),nspec(2),num_sq))

do i=1,num_sq
!call alloc_complex_svd(f_mat(:,:,i), svf(:,i), uf(:,:,i), vtf(:,:,i))
call schmidt_modes(f_mat(:,:,i), svf(:,i), uf(:,:,i), vtf(:,:,i), units(:,i)) 
end do
!call printvectors(uf1, 'uf1')

!close(14)
!close(15)
!close(16)
!close(17)
!close(18)
!close(19)
!close(20)
!close(21)

!call matrixexp

contains 

subroutine matrixexp

!call printvectors( matmul(mat_bs,conjg(transpose(mat_bs))), 'print bs pi/4')
!do i=1, size(f_mat,1)
!    write(15,*) j,i, (real(f_mat(i,j)), j=1, size(f_mat,1))
!end do

f_size=size(f_mat1,1)
sizeofexp=4*f_size
print*, 'dim of exp', sizeofexp
print*, 'transform_no_bs'

jsaamp=1.0_dp
do j=1,1

f_mat1=f_mat1*jsaamp
f_mat2=f_mat2*jsaamp 

!
101 format (2f10.2) 
102 format (16f10.2) 

mat_sq1=make_squeezer(nspace, nspec(1), 1,2,f_mat1)
write(*,101) real(f_mat1)
print*, 'fmat1'
write(*,102) real(mat_sq1) 
print*, 'mat sq 12' 
!
mat_sq2=make_squeezer(nspace, nspec(1), 3,4,f_mat2)


transform_no_bs=matmul(mat_sq1,mat_sq2)
print*, 'nspec',  nspec
g4_nobs=g4(transform_no_bs,nspec(1))

theta_samples=10
allocate(g4norm(theta_samples))
allocate(theta(theta_samples))

theta=0.0_dp
thetaincr=0.0_dp
do i=1, theta_samples
    
    theta(i)=thetaincr
    thetaincr=thetaincr+0.02*pi

    ! modes 2 & 3
    mat_bs=0.0_dp
    call make_bs(nspace,nspec(1),mat_bs,2,3,theta(i))
    !call printvectors(mat_bs, 'beamsplitter')


!transform_no_bs=matmul(mat_sq1,mat_sq2)
!g4_nobs=g4(transform_no_bs,nspec)

    !call printvectors(transform, 'sq1 * sq2')

    !print*, 'now do beam splitter on modes 2&3'

    transform=matmul(mat_bs,transform_no_bs)
    !call printvectors(transform(1:n,1:n), 'sq1*sq2*BS') 
    g4un_norm= g4(transform, nspec(1))
    g4norm(i)=g4un_norm/g4_nobs
    ! normalise 

    print*, 'theta', theta(i), 'g4 no BS', g4_nobs, 'g4 un',g4un_norm, 'g4norm', g4norm(i)
    write(16,*) theta(i), g4norm(i), g4_nobs,g4un_norm
end do ! theta 
    
g4norm=g4norm/(sum(g4norm))

do k=1, size(g4norm)
    write(17,*) jsaamp, theta(k), g4norm(k)
end do

deallocate(g4norm)
deallocate(theta)

! 1 is no scaling
    jsaamp=jsaamp+0.1_dp

end do ! jsa amp

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



end subroutine matrixexp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

end program num_hom
