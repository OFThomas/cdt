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

!subroutine make_sq(nspace,nspec,symp_mat,m1,m2,f_jsa)



!end subroutine make_sq


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
