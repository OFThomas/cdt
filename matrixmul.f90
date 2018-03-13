!Oliver Thomas 2018 Bristol
program matrixmul
implicit none

integer, parameter :: dp=selected_real_kind(15,300)
integer :: n, qubits, numofdecomp, i, j, counter
integer, allocatable, dimension(:) :: p

real(kind=dp), parameter :: r2=1/sqrt(real(2,kind=dp)) 

real(kind=dp), allocatable, dimension(:,:) :: unitary, ident, uprod
real(kind=dp), allocatable, dimension(:,:,:)  :: u

qubits=2
n= 2**qubits
numofdecomp=int(n*(n-1)/2.0_dp)

allocate(ident(n,n))
allocate(unitary(n,n))
allocate(u(n,n,n*n))
allocate(uprod(n,n))
allocate(p(n))

p(1:n)=(/1,2,4,3/)

ident=0.0_dp
unitary=0.0_dp
u=0.0_dp

!make identity
ident=identity(n)

!make u's ident 
do i=1, size(u,3)
  do j=1,n
    u(j,j,i)=1.0_dp
  end do
end do

do i=1, size(u,3)
  u(:,:,i)=identity(n)
end do

!init uprod as ident
uprod=identity(n)


!write(*,*) int(ident)

!make unitary
unitary(1:n,1)=(/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
unitary(1:n,2)=(/0.0_dp, r2, r2, 0.0_dp/)
unitary(1:n,3)=(/0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/)
unitary(1:n,4)=(/0.0_dp, r2, -r2,0.0_dp/)

!write(*,*) unitary

!!!! make unitary gates

!u1
!u(1:n,1,1)=(/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
!u(1:n,2,1)=(/0.0_dp,r2/4.0, 0.0_dp, -r2/4.0/)
!u(1:n,3,1)=(/0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp/)
!u(1:n,4,1)=(/0.0_dp, r2/4.0, 0.0_dp, r2/4.0/)

!u2
!u(1:n,1,2)=(/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
!u(1:n,2,2)=(/0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp/)
!u(1:n,3,2)=(/0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp/)
!u(1:n,4,2)=(/0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp/)

!u3


!print*, p(n), p(n-1), p(1)
!write(*,*) uprod
!write(*,*) u(:,:,1)
uprod=unitary
do i=1,n !col
  do j=1,n-1 !row
    !write(*,*)
    !write(*,*) uprod
	!print*, 
	!print*, p(n+1-j
    if(p(n-j+1).ne.p(i)) then
    call makeunitary(p(n-j),p(n-j+1),p(i), uprod, u(:,:,(j-1)*n+i))
    end if
  end do
end do

print*,
print*, 'unitaries'

do i=1, n*n
  if (icheck(u(:,:,i))==0) then 
print*, i
    write(*,*) u(:,:,i)
    print*,
  end if
end do

print*,
print*, unitarycheck(u,unitary)
if (unitarycheck(u,unitary)==1) then
  print*, 'THIS IS UNITARY'
end if

deallocate(ident)
deallocate(unitary)
deallocate(u)
deallocate(uprod)   
deallocate(p)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

!!!!! Check product gives identity 
function unitarycheck(umatrices, uni)
  real(kind=dp) :: unitarycheck
  real(kind=dp), dimension(:,:,:), intent(in) :: umatrices
  real(kind=dp), dimension(:,:), intent(in) :: uni
  real(kind=dp), dimension(:,:), allocatable :: uprod
  integer :: i, j

unitarycheck=1
uprod=uni
!write(*,*) uprod
!print*,
!do u_n*u_n-1*...*u1*Unitary=Ident
do i=1, size(umatrices,3)-1
  if (icheck(umatrices(:,:,n*n-i))==0) then
    uprod=matmul(uprod(:,:), umatrices(:,:,n*n-i))
    !write(*,*) uprod
    !write(*,*)
  end if
end do

unitarycheck=icheck(uprod)

end function unitarycheck 

!!!! Calc givens rotation
subroutine givensrot(a, b, c, s, r)
  real(kind=dp) :: a, b, c, s, r, h, d
h=0.0
d=0.0
!print*,
 !write(*,*) 'a',a,'b',b,'c',c,'s',s
if (b.ne.0.0_dp) then
  h=hypot(a,b)
  d=1.0_dp/h
  c=abs(a)*d
  s=sign(d,a)*b
  r=sign(1.0_dp,a)*h
else
  c=1.0_dp
  s=0.0_dp
  r=a
  
end if
end subroutine givensrot  

!!!! find type of u
subroutine makeunitary(row1,row2,col,ucurrent, ugate)
  real(kind=dp) :: c,s,r
  real(kind=dp), dimension(:,:), intent(inout) :: ucurrent, ugate
  integer, intent(in) :: row1, row2, col
  
  ugate=identity(size(ucurrent,1))
  c=0.0
  s=0.0
  r=0.0

!print*, 'r1',row1,'r2',row2,'col',col
call givensrot(ucurrent(col,row1), ucurrent(col,row2), c,s,r)
!print*, 'c=', c, 's=', s
  ugate(row1,row1)=c
  ugate(row1,row2)=-s
  ugate(row2,row1)=s
  ugate(row2,row2)=c
 if (icheck(ugate)==0) then 
  print*, 'ugate'
  write(*,*) ugate
  print*,
end if
ucurrent=matmul(ucurrent,ugate)

end subroutine makeunitary

!!!!!! make identiy matrix dim n
function identity(n)
  real(kind=dp), dimension(n,n) :: identity
  integer :: n, i

identity=0.0_dp
do i=1, n
  identity(i,i) =1.0_dp
end do
end function identity


!!!!! Check product gives identity 
function icheck(uni)
  real(kind=dp) :: icheck
  real(kind=dp), dimension(:,:), intent(in) :: uni
  integer :: i, j

icheck=1

!check ident
itest:do i=1, size(uni,1)
  do j=1, size(uni,1)
    if (i.ne.j) then
      if (uni(i,j).ne.0.0_dp) then
        !print*, 'not ident'
        icheck=0
        exit itest
      end if
    else if (i.eq.j) then
      if (uni(i,j).ne.1.0_dp) then
        !print*, 'not identity'
        icheck=0
        exit itest
      end if
    end if
  end do
end do itest

end function icheck



end program matrixmul




