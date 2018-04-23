!Oliver Thomas 2018 Bristol
program matrixmul
implicit none

integer, parameter :: dp=selected_real_kind(15,300)
integer :: n, qubits, numofdecomp, i, j, counter
integer, allocatable, dimension(:) :: p, gatenum

real(kind=dp), parameter :: invr2=1/sqrt(real(2,kind=dp)), invr3=1/sqrt(real(3,kind=dp)), invr6=1/sqrt(real(6,kind=dp))
real(kind=dp), parameter :: r2=sqrt(real(2,kind=dp)) 

real(kind=dp), allocatable, dimension(:,:) :: unitary, ident, uprod
real(kind=dp), allocatable, dimension(:,:,:)  :: u, gateseq

counter=1
print*, 'Enter number of qubits, 2, 3 or 4'
read*, qubits
n= 2**qubits
numofdecomp=int(n*(n-1)/2.0_dp)

allocate(ident(n,n))
allocate(unitary(n,n))
allocate(u(n,n,n*n))
allocate(uprod(n,n))
allocate(p(n))
allocate(gateseq(n,n,n*n))
allocate(gatenum(numofdecomp))


ident=0.0_dp
unitary=0.0_dp
u=0.0_dp
gateseq=0.0_dp
gatenum=1

!#make identity
ident=identity(n)

!#make u's ident 
do i=1, size(u,3)
  u(:,:,i)=identity(n)
end do
gateseq=u
!#init uprod as ident

!#make unitary
if (qubits==2) then 
  p(1:n)=(/1,2,4,3/)

  unitary(1:n,1)=(/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
  unitary(1:n,2)=(/0.0_dp, invr2, invr2, 0.0_dp/)
  unitary(1:n,3)=(/0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/)
  unitary(1:n,4)=(/0.0_dp, invr2, -invr2,0.0_dp/)

else if (qubits==3) then

  p(1:n)=(/1,2,4,3,7,8,6,5/)
!# col,row
  unitary(1,1)=1.0_dp

  unitary(2,2)=invr3
  unitary(2,5)=r2*invr3

  unitary(3,2)=invr3
  unitary(3,5)=-invr6
  unitary(3,7)=invr2

  unitary(4,3)=invr3
  unitary(4,6)=invr6
  unitary(4,8)=invr2

  unitary(5,2)=invr3
  unitary(5,5)=-invr6
  unitary(5,7)=-invr2

  unitary(6,3)=invr3
  unitary(6,6)=invr6
  unitary(6,8)=-invr2

  unitary(7,3)=invr3
  unitary(7,6)=-r2*invr3

  unitary(8,4)=1.0_dp

else if (qubits==4) then

  p(1:n)=(/1,2,4,3,7,8,6,5,13,9,10,12,11,15,16,14/)
  !# col,row
  !0000
  unitary(1,1)=1.0_dp
  !0001
  unitary(2,2)=0.5_dp
  unitary(2,6)=sqrt(0.75_dp)
  !0010
  unitary(3,2)=0.5_dp
  unitary(3,6)=-sqrt(1.0_dp/12.0_dp)
  unitary(3,9)=sqrt(2.0_dp/3.0_dp)
  !0011
  unitary(4,3)=sqrt(1.0_dp/6.0_dp)
  unitary(4,7)=sqrt(1.0_dp/6.0_dp)
  unitary(4,10)=sqrt(1.0_dp/3.0_dp)
  unitary(4,15)=sqrt(1.0_dp/3.0_dp)
  !0100
  unitary(5,2)=0.5_dp
  unitary(5,6)=-sqrt(1.0_dp/12.0_dp)
  unitary(5,9)=-sqrt(1.0_dp/6.0_dp)
  unitary(5,12)=sqrt(1.0_dp/2.0_dp)
  !0101
  unitary(6,3)=sqrt(1.0_dp/6.0_dp)
  unitary(6,7)=sqrt(1.0_dp/6.0_dp)
  unitary(6,10)=-sqrt(1.0_dp/12.0_dp)
  unitary(6,13)=0.5_dp
  unitary(6,15)=-sqrt(1.0_dp/12.0_dp)
  unitary(6,16)=0.5_dp 
  !0110
  unitary(7,3)=sqrt(1.0_dp/6.0_dp)
  unitary(7,7)=-sqrt(1.0_dp/6.0_dp)
  unitary(7,10)=sqrt(1.0_dp/12.0_dp)
  unitary(7,13)=0.5_dp
  unitary(7,15)=-sqrt(1.0_dp/12.0_dp)
  unitary(7,16)=-0.5_dp
  !0111
  unitary(8,4)=0.5_dp
  unitary(8,8)=sqrt(1.0_dp/12.0_dp)
  unitary(8,11)=sqrt(1.0_dp/6.0_dp)
  unitary(8,14)=-sqrt(0.5_dp)
  !1000
  unitary(9,2)=0.5_dp
  unitary(9,6)=-sqrt(1.0_dp/12.0_dp)
  unitary(9,9)=-sqrt(1.0_dp/6.0_dp)
  unitary(9,12)=-sqrt(0.5_dp)
  !1001
  unitary(10,3)=sqrt(1.0_dp/6.0_dp)
  unitary(10,7)=sqrt(1.0_dp/6.0_dp)
  unitary(10,10)=-sqrt(1.0_dp/12.0_dp)
  unitary(10,13)=-0.5_dp
  unitary(10,15)=-sqrt(1.0_dp/12.0_dp)
  unitary(10,16)=-0.5_dp
  !1010
  unitary(11,3)=sqrt(1.0_dp/6.0_dp)
  unitary(11,7)=-sqrt(1.0_dp/6.0_dp)
  unitary(11,10)=sqrt(1.0_dp/12.0_dp)
  unitary(11,13)=-0.5_dp
  unitary(11,15)=-sqrt(1.0_dp/12.0_dp)
  unitary(11,16)=0.5_dp
  !1011
  unitary(12,4)=0.5_dp
  unitary(12,8)=sqrt(1.0_dp/12.0_dp)
  unitary(12,11)=sqrt(1.0_dp/6.0_dp)
  unitary(12,14)=sqrt(0.5)
  !1100
  unitary(13,3)=sqrt(1.0_dp/6.0_dp)
  unitary(13,7)=-sqrt(1.0_dp/6.0_dp)
  unitary(13,10)=-sqrt(1.0_dp/3.0_dp)
  unitary(13,15)=sqrt(1.0_dp/3.0_dp)
  !1101
  unitary(14,4)=0.5_dp
  unitary(14,8)=sqrt(1.0_dp/12.0_dp)
  unitary(14,11)=-sqrt(2.0_dp/3.0_dp)
  !1110
  unitary(15,4)=0.5_dp
  unitary(15,8)=-sqrt(3.0_dp/4.0_dp)
  !1111
  unitary(16,5)=1.0_dp

end if

21 format ( 4F7.3)
write(*,21) unitary
!#!!! make unitary gates

print*,
write(*,21) matmul(unitary,transpose(unitary))

uprod=unitary
do i=1,n !#col
  do j=1,n-1 !#row
    if(p(n-j+1).ne.p(i)) then
      call makeunitary(p(n-j),p(n-j+1),p(i), uprod, u(:,:,(i-1)*n+j))
    end if
  end do
end do

print*,
print*, 'unitaries'

uprod=unitary
do i=1, n*n
  if (icheck(u(:,:,i))==0) then 
    uprod=matmul(uprod(:,:),u(:,:,i))
  end if
end do

  print*, 'THIS IS UNITARY'
  call invert(u,gateseq,counter)
  call gateset(gateseq)
 
  write(*,21) gateseq(:,:,1:counter)
do i=1, counter
  uprod=matmul(uprod,gateseq(:,:,i))
end do

print*, '------------------------------------------------------------------'
print*, 'Unitary matrix from', counter-1, 'gates'
print*, '------------------------------------------------------------------'

print*, size(gateseq,3)

deallocate(ident)
deallocate(unitary)
deallocate(u)
deallocate(uprod)   
deallocate(p)
deallocate(gateseq)
deallocate(gatenum)

!#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

!# print non identity elements
subroutine gateset(matrix)
  real(kind=dp), dimension(:,:,:), intent(in) :: matrix
  integer :: n, i, j, k
  n=size(matrix,3)
  do i=1, n
    if (icheck(matrix(:,:,i))==0) then 
    end if
    if (cnotcheck(matrix(:,:,i))==0) then
    end if
  end do
end subroutine gateset

!# transpose and invert array
subroutine invert(matrix,inverted,count)
  real(kind=dp), dimension(:,:,:), intent(inout) :: inverted
  real(kind=dp), dimension(:,:,:), intent(in) :: matrix
  integer :: i, n, count
  n=size(matrix,3)
  count=1
  do i=1,n
    if (icheck(matrix(:,:,n-i+1))==0) then
      inverted(:,:,count) = transpose(matrix(:,:,n-i+1)) 
      count=count+1
    end if
  end do
  
end subroutine invert

!#!!!! Check product gives identity 
function cnotcheck(uni)
  real(kind=dp) :: cnotcheck
  real(kind=dp), dimension(:,:), intent(in) :: uni
  integer :: i, j

!check for CNOT
cnotcheck=1
!#check if cnot
cnottest:do i=1, size(uni,1)
  do j=1, size(uni,1)
    if ((abs(uni(i,j)).eq.0).or.(abs(uni(i,j)).eq.1)) then
     ! print*, ' 0 or 1'
    else 
      print*, 'not 0 or 1!!!!!'
    end if
  end do
end do cnottest
end function cnotcheck


!#!!!! Check product gives identity 
function icheck(uni)
  real(kind=dp) :: icheck
  real(kind=dp), dimension(:,:), intent(in) :: uni
  integer :: i, j

icheck=1
!#check ident
itest:do i=1, size(uni,1)
  do j=1, size(uni,1)
    if (i.ne.j) then
      if (abs(uni(i,j))>=1e-10) then
        icheck=0
        exit itest
      end if
    else if (i.eq.j) then
      if (abs(abs(uni(i,j))-1)>=1e-10) then
        icheck=0
        exit itest
      end if
    end if
  end do
end do itest
end function icheck

!#!!!! Check product gives identity 
function unitarycheck(umatrices, uni)
  real(kind=dp) :: unitarycheck
  real(kind=dp), dimension(:,:,:), intent(in) :: umatrices
  real(kind=dp), dimension(:,:), intent(in) :: uni
  real(kind=dp), dimension(:,:), allocatable :: uprod
  integer :: i, j

unitarycheck=1
uprod=uni
!#do u_n*u_n-1*...*u1*Unitary=Ident
do i=1, size(umatrices,3)-1
  if (icheck(umatrices(:,:,n*n-i))==0) then
    uprod=matmul(uprod(:,:), umatrices(:,:,i))
  end if
end do

unitarycheck=icheck(uprod)

end function unitarycheck 

!#!!! find type of u
subroutine makeunitary(row1,row2,col,ucurrent, ugate)
  real(kind=dp) :: c,s,r
  real(kind=dp), dimension(:,:), intent(inout) :: ucurrent, ugate
  integer, intent(in) :: row1, row2, col
  
  ugate=identity(size(ucurrent,1))
  c=0.0
  s=0.0
  r=0.0

call givensrot(ucurrent(col,row1), ucurrent(col,row2), c,s,r)
  ugate(row1,row1)=c
  ugate(row1,row2)=-s
  ugate(row2,row1)=s
  ugate(row2,row2)=c
ucurrent=matmul(ucurrent,ugate)
end subroutine makeunitary

!#!!! Calc givens rotation
subroutine givensrot(a, b, c, s, r)
  real(kind=dp) :: a, b, c, s, r, h, d
h=0.0
d=0.0
 !#write(*,*) 'a',a,'b',b,'c',c,'s',s
if (abs(b)>=1e-1) then
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


!#!!!!! make identiy matrix dim n
function identity(n)
  real(kind=dp), dimension(n,n) :: identity
  integer :: n, i

identity=0.0_dp
do i=1, n
  identity(i,i) =1.0_dp
end do
end function identity

end program matrixmul

