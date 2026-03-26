module set_bc_common
  use mod_globals, only : nx, ny, nz
  implicit none
  interface set_bc_cyclic
    module procedure set_bc_cyclic2, set_bc_cyclic4, set_bc_cyclic6
  end interface set_bc_cyclic
contains
  subroutine set_bc_cyclic2(id_accuracy, nx, ny, nz, Q)
    integer(kind=2), intent(in) :: id_accuracy
    integer, intent(in)         :: nx, ny, nz
    real(8), intent(inout)      :: Q(5,nx,ny,nz)
    integer i, j, k, l
    do k = 2, nz-1
      do j = 2, ny-1
        do l = 1, 5
          Q(l,1,j,k)  = Q(l,nx-1,j,k)
          Q(l,nx,j,k) = Q(l,2,j,k)
    enddo;enddo;enddo
    do k = 2, nz-1
      do i = 2, nx-1
        do l = 1, 5
          Q(l,i,1,k)  = Q(l,i,ny-1,k)
          Q(l,i,ny,k) = Q(l,i,2,k)
    enddo;enddo;enddo
    do k = 2, nz-1
      do l = 1, 5
        Q(l,1,1,k)   = Q(l,nx-1,ny-1,k)
        Q(l,nx,1,k)  = Q(l,2,ny-1,k)
        Q(l,1,ny,k)  = Q(l,nx-1,2,k)
        Q(l,nx,ny,k) = Q(l,2,2,k)
    enddo;enddo
    do j = 1, ny
      do i = 1, nx
        do l = 1, 5
          Q(l,i,j,1)  = Q(l,i,j,nz-1)
          Q(l,i,j,nz) = Q(l,i,j,2)
    enddo;enddo;enddo
  end subroutine set_bc_cyclic2


  subroutine set_bc_cyclic4(id_accuracy, nx, ny, nz, Q)
    integer(kind=4), intent(in) :: id_accuracy
    integer, intent(in)         :: nx, ny, nz
    real(8), intent(inout)      :: Q(5,nx,ny,nz)
    integer i, j, k, l
    do k = 3, nz-2
      do j = 3, ny-2
        do l = 1, 5
          Q(l,1,j,k) = Q(l,nx-3,j,k)
          Q(l,2,j,k) = Q(l,nx-2,j,k)
          Q(l,nx-1,j,k) = Q(l,3,j,k)
          Q(l,nx,j,k) = Q(l,4,j,k)
    enddo;enddo;enddo
    do k = 3, nz-2
      do i = 3, nx-2
        do l = 1, 5
          Q(l,i,1,k) = Q(l,i,ny-3,k)
          Q(l,i,2,k) = Q(l,i,ny-2,k)
          Q(l,i,ny-1,k) = Q(l,i,3,k)
          Q(l,i,ny,k) = Q(l,i,4,k)
    enddo;enddo;enddo
    do k = 3, nz-2
      do l = 1, 5
        Q(l,1,1,k) = Q(l,nx-3,ny-3,k)
        Q(l,1,2,k) = Q(l,nx-3,ny-2,k)
        Q(l,2,1,k) = Q(l,nx-2,ny-3,k)
        Q(l,2,2,k) = Q(l,nx-2,ny-2,k)
        Q(l,nx-1,1,k) = Q(l,3,ny-3,k)
        Q(l,nx-1,2,k) = Q(l,3,ny-2,k)
        Q(l,nx,1,k)   = Q(l,4,ny-3,k)
        Q(l,nx,2,k)   = Q(l,4,ny-2,k)
        Q(l,1,ny-1,k) = Q(l,nx-3,3,k)
        Q(l,1,ny,k)   = Q(l,nx-3,4,k)
        Q(l,2,ny-1,k) = Q(l,nx-2,3,k)
        Q(l,2,ny,k)   = Q(l,nx-2,4,k)
        Q(l,nx-1,ny-1,k) = Q(l,3,3,k)
        Q(l,nx-1,ny,k)   = Q(l,3,4,k)
        Q(l,nx,ny-1,k)   = Q(l,4,3,k)
        Q(l,nx,ny,k)     = Q(l,4,4,k)
    enddo;enddo
    do j = 1, ny
      do i = 1, nx
        do l = 1, 5
          Q(l,i,j,1) = Q(l,i,j,nz-3)
          Q(l,i,j,2) = Q(l,i,j,nz-2)
          Q(l,i,j,nz-1) = Q(l,i,j,3)
          Q(l,i,j,nz) = Q(l,i,j,4)
    enddo;enddo;enddo
  end subroutine set_bc_cyclic4


  subroutine set_bc_cyclic6(id_accuracy, nx, ny, nz, Q)
    integer(kind=8), intent(in) :: id_accuracy
    integer, intent(in)         :: nx, ny, nz
    real(8), intent(inout)      :: Q(5,nx,ny,nz)
    integer i, j, k, l
    do k = 4, nz-3
      do j = 4, ny-3
        do l = 1, 5
          Q(l,1,j,k) = Q(l,nx-5,j,k)
          Q(l,2,j,k) = Q(l,nx-4,j,k)
          Q(l,3,j,k) = Q(l,nx-3,j,k)
          Q(l,nx-2,j,k) = Q(l,4,j,k)
          Q(l,nx-1,j,k) = Q(l,5,j,k)
          Q(l,nx,j,k)   = Q(l,6,j,k)
    enddo;enddo;enddo
    do k = 4, nz-3
      do i = 4, nx-3
        do l = 1, 5
          Q(l,i,1,k) = Q(l,i,ny-5,k)
          Q(l,i,2,k) = Q(l,i,ny-4,k)
          Q(l,i,3,k) = Q(l,i,ny-3,k)
          Q(l,i,ny-2,k) = Q(l,i,4,k)
          Q(l,i,ny-1,k) = Q(l,i,5,k)
          Q(l,i,ny,k)   = Q(l,i,6,k)
    enddo;enddo;enddo
    do k = 4, nz-3
      do l = 1, 5
        Q(l,1,1,k) = Q(l,nx-5,ny-5,k)
        Q(l,1,2,k) = Q(l,nx-5,ny-4,k)
        Q(l,1,3,k) = Q(l,nx-5,ny-3,k)
        Q(l,2,1,k) = Q(l,nx-4,ny-5,k)
        Q(l,2,2,k) = Q(l,nx-4,ny-4,k)
        Q(l,2,3,k) = Q(l,nx-4,ny-3,k)
        Q(l,3,1,k) = Q(l,nx-3,ny-5,k)
        Q(l,3,2,k) = Q(l,nx-3,ny-4,k)
        Q(l,3,3,k) = Q(l,nx-3,ny-3,k)
        Q(l,nx-2,1,k) = Q(l,4,ny-5,k)
        Q(l,nx-2,2,k) = Q(l,4,ny-4,k)
        Q(l,nx-2,3,k) = Q(l,4,ny-3,k)
        Q(l,nx-1,1,k) = Q(l,5,ny-5,k)
        Q(l,nx-1,2,k) = Q(l,5,ny-4,k)
        Q(l,nx-1,3,k) = Q(l,5,ny-3,k)
        Q(l,nx,1,k)   = Q(l,6,ny-5,k)
        Q(l,nx,2,k)   = Q(l,6,ny-4,k)
        Q(l,nx,3,k)   = Q(l,6,ny-3,k)
        Q(l,1,ny-2,k) = Q(l,nx-5,4,k)
        Q(l,1,ny-1,k) = Q(l,nx-5,5,k)
        Q(l,1,ny,k)   = Q(l,nx-5,6,k)
        Q(l,2,ny-2,k) = Q(l,nx-4,4,k)
        Q(l,2,ny-1,k) = Q(l,nx-4,5,k)
        Q(l,2,ny,k)   = Q(l,nx-4,6,k)
        Q(l,3,ny-2,k) = Q(l,nx-3,4,k)
        Q(l,3,ny-1,k) = Q(l,nx-3,5,k)
        Q(l,3,ny,k)   = Q(l,nx-3,6,k)
        Q(l,nx-2,ny-2,k) = Q(l,4,4,k)
        Q(l,nx-2,ny-1,k) = Q(l,4,5,k)
        Q(l,nx-2,ny,k)   = Q(l,4,6,k)
        Q(l,nx-1,ny-2,k) = Q(l,5,4,k)
        Q(l,nx-1,ny-1,k) = Q(l,5,5,k)
        Q(l,nx-1,ny,k)   = Q(l,5,6,k)
        Q(l,nx,ny-2,k)   = Q(l,6,4,k)
        Q(l,nx,ny-1,k)   = Q(l,6,5,k)
        Q(l,nx,ny,k)     = Q(l,6,6,k)
    enddo;enddo
    do j = 1, ny
      do i = 1, nx
        do l = 1, 5
          Q(l,i,j,1) = Q(l,i,j,nz-5)
          Q(l,i,j,2) = Q(l,i,j,nz-4)
          Q(l,i,j,3) = Q(l,i,j,nz-3)
          Q(l,i,j,nz-2) = Q(l,i,j,4)
          Q(l,i,j,nz-1) = Q(l,i,j,5)
          Q(l,i,j,nz)   = Q(l,i,j,6)
    enddo;enddo;enddo
  end subroutine set_bc_cyclic6
end module set_bc_common
