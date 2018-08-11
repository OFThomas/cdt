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

subroutine make_sq(nspace,nspec,symp_mat,m1,m2, alpha, beta)

    integer :: nspace, nspec, m1, m2
    complex(kind=dp), dimension(:,:), allocatable :: symp_mat, alpha, beta
    
    !counters
    integer :: i,j, n
    integer :: m1s, m2s

    n=nspace*nspec
    m1s=m1*nspec-1
    m2s=m2*nspec-1

    !if ((m1s <= i < m1s + nspec) .or. (m2s <= i < m2s + nspec)) then
    !>@note alpha & beta are 2 spatial modes and all spectral modes
    !> dim 2*nspace*nspec
   
    ! make passives identity 
    symp_mat(:,:)=c_identity(n*2)

    !then overwrite values from alpha
    ! do alpha first on modes m1, m2 
    ! for m1 do for all spectral modes copy alpha
    symp_mat(m1s:m1s+nspec, m1s:m1s+nspec)=alpha(1:nspec, 1:nspec)
    ! for m2 do the same
    symp_mat(m2s:m2s+nspec, m2s:m2s+nspec)=alpha(nspec+1:2*nspec, nspec+1:2*nspec)

    ! alpha conjg is then 
    symp_mat(n+1:2*n, n+1:2*n)=conjg(symp_mat(1:n, 1:n))

    ! beta is much easier
    !same as alpha but second index is off-set by n
   

    ! for mode 1
    !>@TODO check this is legal...
    !> full diag sq
    !> symp_mat(m1s:m1s+nspec, m1s+n:m1s+nspec+n)=beta(1:nspec, 1+nspec:2*nspec)
   
    ! mode 1
    symp_mat(m1s:m1s+nspec, m2s+n:m2s+nspec+n)=beta(1:nspec, 1+nspec:2*nspec)
    
    ! for mode 2
    symp_mat(m2s:m2s+nspec, m1s+n:m1s+nspec+n)=beta(nspec+1:2*nspec, 1:nspec)
    
    !>@TODO
    !> probably not legal
    !>symp_mat(m2s:m2s+nspec, m2s+n:m2s+nspec+n)=beta(nspec+1:2*nspec, 1:nspec)

    ! beta conjg is then
    symp_mat(n+1:2*n, 1:n)=conjg(symp_mat(1:n, n+1:2*n))


end subroutine make_sq


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




end module makeopticalelements
