module set_bc_common
  implicit none
contains
  subroutine set_bc_cyclic(nx, ny, nz, Q)
    integer, intent(in), value :: nx, ny, nz
    real(8), intent(inout)     :: Q(5,nx,ny,nz)
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
  end subroutine set_bc_cyclic
end module set_bc_common
