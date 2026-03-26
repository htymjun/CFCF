module set_bc_common
  use mod_globals, only : nx, ny, nz
  implicit none
  interface set_bc_cyclic
    module procedure set_bc_cyclic2_init, set_bc_cyclic4_init, set_bc_cyclic6_init
  end interface set_bc_cyclic
contains
  subroutine set_bc_cyclic2_init(id_accuracy, nx, ny, nz, Q)
    integer(kind=2), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(inout)             :: Q(5,nx,ny,nz)
    integer i, j, k
    do k = 2, nz-1
      do j = 2, ny-1
        Q(:,1,j,k)  = Q(:,nx-1,j,k)
        Q(:,nx,j,k) = Q(:,2,j,k)
    enddo;enddo
    do k = 2, nz-1
      do i = 2, nx-1
        Q(:,i,1,k)  = Q(:,i,ny-1,k)
        Q(:,i,ny,k) = Q(:,i,2,k)
    enddo;enddo
    do k = 2, nz-1
      Q(:,1,1,k)   = Q(:,nx-1,ny-1,k)
      Q(:,nx,1,k)  = Q(:,2,ny-1,k)
      Q(:,1,ny,k)  = Q(:,nx-1,2,k)
      Q(:,nx,ny,k) = Q(:,2,2,k)
    enddo
    do j = 1, ny
      do i = 1, nx
        Q(:,i,j,1)  = Q(:,i,j,nz-1)
        Q(:,i,j,nz) = Q(:,i,j,2)
    enddo;enddo
  end subroutine set_bc_cyclic2_init 


  subroutine set_bc_cyclic4_init(id_accuracy, nx, ny, nz, Q)
    integer(kind=4), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(inout)             :: Q(5,nx,ny,nz)
    integer i, j, k
    do k = 3, nz-2
      do j = 3, ny-2
        Q(:,1:2,j,k) = Q(:,nx-3:nx-2,j,k)
        Q(:,nx-1:nx,j,k) = Q(:,3:4,j,k)
    enddo;enddo
    do k = 3, nz-2
      do i = 3, nx-2
        Q(:,i,1:2,k) = Q(:,i,ny-3:ny-2,k)
        Q(:,i,ny-1:ny,k) = Q(:,i,3:4,k)
    enddo;enddo
    do k = 3, nz-2
      Q(:,1:2,1:2,k) = Q(:,nx-3:nx-2,ny-3:ny-2,k)
      Q(:,nx-1:nx,1:2,k) = Q(:,3:4,ny-3:ny-2,k)
      Q(:,1:2,ny-1:ny,k) = Q(:,nx-3:nx-2,3:4,k)
      Q(:,nx-1:nx,ny-1:ny,k) = Q(:,3:4,3:4,k)
    enddo
    do j = 1, ny
      do i = 1, nx
        Q(:,i,j,1:2) = Q(:,i,j,nz-3:nz-2)
        Q(:,i,j,nz-1:nz) = Q(:,i,j,3:4)
    enddo;enddo
  end subroutine set_bc_cyclic4_init


  subroutine set_bc_cyclic6_init(id_accuracy, nx, ny, nz, Q)
    integer(kind=8), intent(in), value :: id_accuracy
    integer, intent(in), value         :: nx, ny, nz
    real(8), intent(inout)             :: Q(5,nx,ny,nz)
    integer i, j, k
    do k = 4, nz-3
      do j = 4, ny-3
        Q(:,1:3,j,k) = Q(:,nx-5:nx-3,j,k)
        Q(:,nx-2:nx,j,k) = Q(:,4:6,j,k)
    enddo;enddo
    do k = 4, nz-3
      do i = 4, nx-3
        Q(:,i,1:3,k) = Q(:,i,ny-5:ny-3,k)
        Q(:,i,ny-2:ny,k) = Q(:,i,4:6,k)
    enddo;enddo
    do k = 4, nz-3
      Q(:,1:3,1:3,k) = Q(:,nx-5:nx-3,ny-5:ny-3,k)
      Q(:,nx-2:nx,1:3,k) = Q(:,4:6,ny-5:ny-3,k)
      Q(:,1:3,ny-2:ny,k) = Q(:,nx-5:nx-3,4:6,k)
      Q(:,nx-2:nx,ny-2:ny,k) = Q(:,4:6,4:6,k)
    enddo
    do j = 1, ny
      do i = 1, nx
        Q(:,i,j,1:3) = Q(:,i,j,nz-5:nz-3)
        Q(:,i,j,nz-2:nz) = Q(:,i,j,4:6)
    enddo;enddo
  end subroutine set_bc_cyclic6_init



end module set_bc_common

