!>@brief program to compute matrix of a JSA
program num_hom
use olis_f90stdlib
implicit none

!>@param dp 
integer, parameter :: dp=selected_real_kind(15,300)

!>@param i counter
!>@param j counter
!>@param k counter
integer :: i,j,k

integer :: w1_steps, w2_steps

!>@param w1 signal freq
!>@param w2 idler freq
!>@param sigma is gaussian width
real(kind=dp) :: w1, w2, sigma

! loop params
real(kind=dp) :: w1_start, w1_end, w2_start, w2_end
real(kind=dp) :: w1_incr, w2_incr
!>@param f_mat matrix for values of function, f
complex(kind=dp), dimension(:,:), allocatable :: f_mat

open(unit=15,file='fplot.dat', status='replace')

w1_start=-2.0_dp
w2_start=-2.0_dp

w1_end=-w1_start
w2_end=-w2_start

w1_incr=0.05_dp
w2_incr=0.05_dp

w1_steps=ceiling((w1_end-w1_start)/w1_incr)
w2_steps=ceiling((w2_end-w2_start)/w2_incr)

sigma=1.0_dp

allocate(f_mat(w1_steps,w2_steps))

w2=w2_start
do j=1,w2_steps
    w1=w1_start
    do i= 1,w1_steps 
        write(15,*) w1, w2, real(f(w1,w2,sigma),kind=dp)
        f_mat(i,j)=real(f(w1,w2,sigma),kind=dp)
        w1=w1+w1_incr
    end do
    w2=w2+w2_incr
end do

close(15)
contains 

!>@brief JSA function taking two freq
!>@param w1 input signal freq
!>@param w2 input idler freq
!>@param sig input variance
    function f(w1,w2, sig)
    complex(kind=dp) :: f
    real(kind=dp), intent(in) :: w1,w2, sig
    
    f=(1.0_dp/(sig*sqrt(2.0_dp*pi)))**2 * exp(-0.5_dp*(w1/sig)**2)*exp(-0.5_dp*(w2/sig)**2)

    end function f
end program num_hom
