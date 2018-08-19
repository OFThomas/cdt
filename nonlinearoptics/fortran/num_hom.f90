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
real(kind=dp) :: w1, w2, sigma1, sigma2

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
real(kind=dp) :: thetaincr, g4_nobs, g4un_norm, jsaamp
real(kind=dp), dimension(:), allocatable :: g4norm, theta

integer :: theta_samples


!!!!
complex(kind=dp), dimension(:,:), allocatable :: f_mat1, f_mat2

!>@param uf1 is u matrix from f_mat1 svd
!>@param vtf1 is vt matrix from f_mat1 svd
complex(kind=dp), dimension(:,:), allocatable :: uf1, vtf1, uf2, vtf2
!>@param svf1 singular values for f_mat1
real(kind=dp), dimension(:), allocatable :: svf1, svf2

!>@param signalfreq array of schmidt decomp values (u matrix)
real(kind=dp) :: signalfreq
!>@param idlerfreq array of schmidt decomp values (vt matrix)
real(kind=dp) :: idlerfreq

!>@note files to write to  
open(unit=14,file='fplotw1w2.dat', status='replace')
open(unit=15,file='fplotw3w4.dat', status='replace')
open(unit=16, file='g4f90data.dat', status='replace')
open(unit=17, file='g4splot.dat', status='replace')
!>
!> jsa 1 
open(unit=20, file='signalfreq1.dat', status='replace')
open(unit=21, file='idlerfreq1.dat', status='replace')
!> jsa 2
open(unit=22, file='signalfreq2.dat', status='replace')
open(unit=23, file='idlerfreq2.dat', status='replace')

w1_start=-6.0_dp
w2_start=-6.0_dp

w1_end=-w1_start
w2_end=-w2_start

!0.05 
w1_incr=0.10_dp
w2_incr=0.10_dp

sigma1=1.0_dp
sigma2=0.5_dp*sigma1

f_mat1= gen_jsa(w1_start, w1_end, w1_incr, w2_start, w2_end, w2_incr, sigma1, sigma2,14)
f_mat2=gen_jsa(w1_start, w1_end, w1_incr, w2_start, w2_end, w2_incr, sigma1, sigma2,15, w1offset=2.0_dp, w2offset=2.0_dp)

! normalise?
f_mat1=f_mat1/(sum(f_mat1)*w1_incr*w2_incr)
f_mat2=f_mat2/(sum(f_mat2)*w1_incr*w2_incr)
write(*,*) 'sum f1',  sum(f_mat1)*w1_incr*w2_incr,  'sum f2', sum(f_mat2)*w1_incr*w2_incr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

nspace=4
nspec=size(f_mat1,1)
print*, 'nspectral dof', nspec
n=nspace*nspec
allocate(mat_bs(2*n,2*n))

call alloc_temparrays(nspace,nspec)

! jsa to be decomp, singular vals,
! u, vt
!>@brief allocates the singular values, u and vt matrices for svd
print*, 's1 f', size(f_mat1,1), 's2 f', size(f_mat1,2)

write(*,*) size(f_mat1)

call alloc_complex_svd(f_mat1, svf1, uf1, vtf1)

call alloc_complex_svd(f_mat2, svf2, uf2, vtf2)

call complex_svd(f_mat1, svf1, uf1, vtf1)

call complex_svd(f_mat2, svf2, uf2, vtf2)
print*, 'done svd'

103 format (3f10.2) 
write(*,103) svf1
print*, 'singular vals'
!write(*,103) real(f_mat1)
print*, 'u'
!write(*,103) real(uf1)
print*, 'vt'
!write(*,103) real(vtf1)

print*, 'matmul'
!write(*,103) real
!call printvectors(matmul(uf1,matmul(f_mat1,vtf1)))

print*, 'so f = singular*u*vt?'
!write(*,103) real
!call printvectors(matmul(f_mat1,matmul(uf1,vtf1)))


print*, 'find element' 
!>@note returns the w1,w2 element from the Jsa
print*, find_element(3,2,uf1,vtf1, svf1)

print*, 'end of prog?'

call write_sigidler(svf1, uf1, vtf1, 20, 21)

call write_sigidler(svf2, uf2, vtf2, 22, 23)


!close(14)
!close(15)
!close(16)
!close(17)
!close(18)
!close(19)
!close(20)
!close(21)

!>@note after doing svd
!> Unitary  = exp(SUM_k r_k * A^H_k * B^H_k -h.c.)
!>          = X_k exp(r_k * A^H_k * B^H_K -h.c.)
!>          = X_k S^ab_k(-r_k)
!>
!> A_k -> cosh(r_k)A_k + sinh(r_k)B^H_k
!> B_k -> cosh(r_k)B_k + sinh(r_k)A^H_k


!call matrixexp

contains 

! make a list of w1 & signalfreq from schmidt decomp
!>@note k is the k modes from schmidt decomp
!> l is the frequency range 

! only print the non-zero k-th singular vals
subroutine write_sigidler(sv,u,vt, sigout, idlerout)
    real(kind=dp), dimension(:) :: sv
    complex(kind=dp), dimension(:,:) :: u, vt
    integer :: sigout, idlerout

do k=1, size(u,1)
    if (abs(sv(k)) >= 1e-10) then
    do l=1, size(u,2)
        !print*, 'k', k, 'l', l
        write(sigout,*)k,l, calc_sig(k,l,u)
    end do
! make a list of w2 and idlerfreq from schmidt decomp
!>@note k is the k modes from schmidt decomp
!> l is the frequency range 
    do l=1, size(vt,1) 
        !print*, 'k', k, 'l', l
        write(idlerout,*) k,l, calc_idler(k,l,vt)
    end do
end if
end do
end subroutine write_sigidler

!>@brief does the SVD of the Jsa 
!>@detail takes w1, w2 and u, vt, sv from SVD and returns the f(w1,w2) element
!> A * f(w1,w2) = SUM_k (r_k*Psi_k(w1)*Phi_k(w2)
!> the k-th row and w1-th column of PSI
!> the w2-th row and k-th column of PHI
!>
!> A_k = INT  dw1 * Psi_k(w1)*a_1(w1)
!> which is integral u(w1,k) 
!>
!> B_k = INT dw2 * Phis_k(w2)*a_2(w2)
!> which is integral vt(k,w2) 
function find_element(w1,w2,a,b, sv)
complex(kind=dp) :: find_element
integer, intent(in) :: w1, w2
integer :: k
complex(kind=dp), dimension(:,:), allocatable, intent(in) :: a, b
real(kind=dp), dimension(:), allocatable, intent(in) :: sv

complex(kind=dp) :: summation 
summation=0.0_dp
!print*, 'test'!size(a,1)
do k=1, size(a,1)
    summation=summation + sv(k) * a(w1,k) * b(k,w2)
!    !print*, 'k', k, 'find el', summation
end do
find_element=summation
end function find_element

function calc_sig(k, ws,u)
real(kind=dp) :: calc_sig
complex(kind=dp), dimension(:,:), intent(in) :: u
integer, intent(in) :: k,ws
calc_sig=u(ws,k)
end function calc_sig

function calc_idler(k, wi, vt)
real(kind=dp) :: calc_idler
complex(kind=dp), dimension(:,:), intent(in) :: vt
integer, intent(in) ::  k, wi
calc_idler=vt(k,wi)
end function calc_idler    

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

mat_sq1=make_squeezer(nspace, nspec, 1,2,f_mat1)
write(*,101) real(f_mat1)
print*, 'fmat1'
write(*,102) real(mat_sq1) 
print*, 'mat sq 12' 
!
mat_sq2=make_squeezer(nspace, nspec, 3,4,f_mat2)


transform_no_bs=matmul(mat_sq1,mat_sq2)
print*, 'nspec',  nspec
g4_nobs=g4(transform_no_bs,nspec)

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
    call make_bs(nspace,nspec,mat_bs,2,3,theta(i))
    !call printvectors(mat_bs, 'beamsplitter')


!transform_no_bs=matmul(mat_sq1,mat_sq2)
!g4_nobs=g4(transform_no_bs,nspec)

    !call printvectors(transform, 'sq1 * sq2')

    !print*, 'now do beam splitter on modes 2&3'

    transform=matmul(mat_bs,transform_no_bs)
    !call printvectors(transform(1:n,1:n), 'sq1*sq2*BS') 
    g4un_norm= g4(transform, nspec)
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
