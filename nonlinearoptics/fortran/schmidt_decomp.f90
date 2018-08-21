!>@brief program to calculate occupied Schmidt-modes of a JSA
module schmidt_decomp
use olis_f90stdlib
use makeopticalelements  
implicit none

private 
!>@param dp 
integer, parameter, private :: dp=selected_real_kind(15,300)

!>@param i counter
!>@param j counter
!>@param k counter
integer :: i,j,k, l

public :: schmidt_modes
contains 

subroutine schmidt_modes(f_mat, svf, uf, vtf, writeout)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>@param uf1 is u matrix from f_mat1 svd
!>@param vtf1 is vt matrix from f_mat1 svd
complex(kind=dp), dimension(:,:) :: uf, vtf
!>@param svf1 singular values for f_mat1
real(kind=dp), dimension(:) :: svf

integer, dimension(:) :: writeout 

complex(kind=dp), dimension(:,:) :: f_mat

!>@note files to write to  
!open(unit=14,file='fplotw1w2.dat', status='replace')
!open(unit=15,file='fplotw3w4.dat', status='replace')
!>
!> jsa 1 
!open(unit=20, file='signalfreq1.dat', status='replace')
!open(unit=21, file='idlerfreq1.dat', status='replace')
!> jsa 2
!open(unit=22, file='signalfreq2.dat', status='replace')
!open(unit=23, file='idlerfreq2.dat', status='replace')

call complex_svd(f_mat, svf, uf, vtf)

print*, 'done svd'

write(*,*) maxval(svf)

print*, 'find element' 
!>@note returns the w1,w2 element from the Jsa
!print*, find_element(3,2,uf1,vtf1, svf1)

!>@note after doing svd
!> Unitary  = exp(SUM_k r_k * A^H_k * B^H_k -h.c.)
!>          = X_k exp(r_k * A^H_k * B^H_K -h.c.)
!>          = X_k S^ab_k(-r_k)
!>
!> A_k -> cosh(r_k)A_k + sinh(r_k)B^H_k
!> B_k -> cosh(r_k)B_k + sinh(r_k)A^H_k


call write_sigidler(svf, uf, vtf, writeout(1), writeout(2))

end subroutine schmidt_modes

! make a list of w1 & signalfreq from schmidt decomp
!>@note k is the k modes from schmidt decomp
!> l is the frequency range 

! only print the non-zero k-th singular vals
subroutine write_sigidler(sv,u,vt, sigout, idlerout)
    real(kind=dp), dimension(:) :: sv
    complex(kind=dp), dimension(:,:) :: u, vt
    integer :: sigout, idlerout

do k=1, size(u,1)
    if (abs(sv(k)) >= 1e-1) then
    print*, sv(k)
        do l=1, size(u,2)
        !print*, 'k', k, 'l', l
        write(sigout,*)k,l, abs(calc_sig(k,l,u))
 
    
    ! make a list of w2 and idlerfreq from schmidt decomp
!>@note k is the k modes from schmidt decomp
!> l is the frequency range 
    !do l=1, size(vt,1) 
        !print*, 'k', k, 'l', l
        write(idlerout,*) k,l, abs(calc_idler(k,l,vt))
    end do
    write(sigout,*) char(10), char(10)
    write(idlerout,*) char(10), char(10)
end if
end do
end subroutine write_sigidler

!>@brief does the SVD of the Jsa 
!>@detail takes w1, w2 and u, vt, sv from SVD and returns the f(w1,w2) element
!> A * f(w1,w2) = SUM_k (r_k*Psi_k(w1)*Phi_k(w2)
!> the w1-th row and k-th column of PSI
!> the k-th row and w2-th column of PHI
!>
!> SUM_k u(w1,k) * vt(k,w2)
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

end module schmidt_decomp
