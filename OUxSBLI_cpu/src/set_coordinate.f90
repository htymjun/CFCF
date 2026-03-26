module set_coordinate
  implicit none
  interface set_grid_cyclic
    module procedure set_grid_cyclic2_3D, set_grid_cyclic4_3D, set_grid_cyclic6_3D
  end interface set_grid_cyclic
contains
  subroutine set_xix(nx, dx, xix)
    integer, intent(in)  :: nx
    real(8), intent(in)  :: dx(nx-1)
    real(8), intent(out) :: xix(nx-1)
    integer i
    do i = 1, nx-1
      xix(i) = 1.d0 / dx(i)
    enddo
  end subroutine set_xix


  subroutine set_etay(ny, dy, etay)
    integer, intent(in)  :: ny
    real(8), intent(in)  :: dy(ny-1)
    real(8), intent(out) :: etay(ny-1)
    integer j
    do j = 1, ny-1
      etay(j) = 1.d0 / dy(j)
    enddo
  end subroutine set_etay


  subroutine set_zetaz(nz, dz, zetaz)
    integer, intent(in)  :: nz
    real(8), intent(in)  :: dz(nz-1)
    real(8), intent(out) :: zetaz(nz-1)
    integer k
    do k = 1, nz-1
      zetaz(k) = 1.d0 / dz(k)
    enddo
  end subroutine set_zetaz


  subroutine set_Jacobian_xy(nx, ny, nz, dx, dy, dz, Jacobian)
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(in)  :: dx(nx-1), dy(ny-1), dz(nz-1)
    real(8), intent(out) :: Jacobian(nx,ny)
    integer i, j, k
    do k = 2, nz-1
      do j = 2, ny-1
        do i = 2, nx-1
          Jacobian(i,j) = 8.d0 / &
          & ((dx(i-1) + dx(i)) * (dy(j-1) + dy(j)) * (dz(k-1) + dz(k)))
    enddo;enddo;enddo
    do i = 1, nx
      Jacobian(i,1)  = Jacobian(i,2)
      Jacobian(i,ny) = Jacobian(i,ny-1)
    enddo
    do j = 1, ny
      Jacobian(1,j)  = Jacobian(2,j)
      Jacobian(nx,j) = Jacobian(nx-1,j)
    enddo
  end subroutine set_Jacobian_xy


  subroutine set_grid_cyclic2_3D(id_accuracy, nx, ny, nz, Lx, Ly, Lz, xc, yc, zc, dx, dy, dz)
    integer(kind=2), intent(in) :: id_accuracy
    integer, intent(in)         :: nx, ny, nz
    real(8), intent(in)         :: Lx, Ly, Lz
    real(8), intent(out)        :: xc(nx), yc(ny), zc(nz), dx(nx-1), dy(ny-1), dz(nz-1)
    real(8) dx1, dy1, dz1, x(nx+1), y(ny+1), z(nz+1)
    integer i, j, k
    x = 0.d0; y = 0.d0; z = 0.d0  ! Initialize arrays
    dx1 = Lx / dble(nx-2)
    dy1 = Ly / dble(ny-2)
    dz1 = Lz / dble(nz-2)
    dx(:) = dx1
    dy(:) = dy1
    dz(:) = dz1
    ! x direction
    do i = 2, nx
      x(i) = dx1 * dble(i-2)
    enddo
    x(1)    = x(2)  - dx1
    x(nx+1) = x(nx) + dx1
    ! y direction
    do j = 2, ny
      y(j) = dy1 * dble(j-2)
    enddo
    y(1)    = y(2)  - dy1
    y(ny+1) = y(ny) + dy1
    ! z direction
    do k = 2, nz
      z(k) = dz1 * dble(k-2)
    enddo
    z(1)    = z(2)  - dz1
    z(nz+1) = z(nz) + dz1
    ! cell centered
    do i = 1, nx
      xc(i) = 0.5d0 * (x(i) + x(i+1))
    enddo
    do j = 1, ny
      yc(j) = 0.5d0 * (y(j) + y(j+1))
    enddo
    do k = 1, nz
      zc(k) = 0.5d0 * (z(k) + z(k+1))
    enddo
  end subroutine set_grid_cyclic2_3D


  subroutine set_grid_cyclic4_3D(id_accuracy, nx, ny, nz, Lx, Ly, Lz, xc, yc, zc, dx, dy, dz)
    integer(kind=4), intent(in) :: id_accuracy
    integer, intent(in)         :: nx, ny, nz
    real(8), intent(in)         :: Lx, Ly, Lz
    real(8), intent(out)        :: xc(nx), yc(ny), zc(nz), dx(nx-1), dy(ny-1), dz(nz-1)
    real(8) dx1, dy1, dz1, x(nx+1), y(ny+1), z(nz+1)
    integer i, j, k
    x = 0.d0; y = 0.d0; z = 0.d0  ! Initialize arrays
    dx1 = Lx / dble(nx-4)
    dy1 = Ly / dble(ny-4)
    dz1 = Lz / dble(nz-4)
    dx(:) = dx1
    dy(:) = dy1
    dz(:) = dz1
    ! x direction
    do i = 3, nx-1
      x(i) = dx1 * dble(i-3)
    enddo
    x(1)    = x(3)    - 2.d0 * dx1
    x(2)    = x(3)    - dx1
    x(nx)   = x(nx-1) + dx1
    x(nx+1) = x(nx-1) + 2.d0 * dx1
    ! y direction
    do j = 3, ny-1
      y(j) = dy1 * dble(j-3)
    enddo
    y(1)    = y(3)    - 2.d0 * dy1
    y(2)    = y(3)    - dy1
    y(ny)   = y(ny-1) + dy1
    y(ny+1) = y(ny-1) + 2.d0 * dy1
    ! z direction
    do k = 3, nz-1
      z(k) = dz1 * dble(k-3)
    enddo
    z(1)    = z(3)    - 2.d0 * dz1
    z(2)    = z(3)    - dz1
    z(nz)   = z(nz-1) + dz1
    z(nz+1) = z(nz-1) + 2.d0 * dz1
    ! cell centered
    do i = 1, nx
      xc(i) = 0.5d0 * (x(i) + x(i+1))
    enddo
    do j = 1, ny
      yc(j) = 0.5d0 * (y(j) + y(j+1))
    enddo
    do k = 1, nz
      zc(k) = 0.5d0 * (z(k) + z(k+1))
    enddo
  end subroutine set_grid_cyclic4_3D


  subroutine set_grid_cyclic6_3D(id_accuracy, nx, ny, nz, Lx, Ly, Lz, xc, yc, zc, dx, dy, dz)
    integer(kind=8), intent(in) :: id_accuracy
    integer, intent(in)         :: nx, ny, nz
    real(8), intent(in)         :: Lx, Ly, Lz
    real(8), intent(out)        :: xc(nx), yc(ny), zc(nz), dx(nx-1), dy(ny-1), dz(nz-1)
    real(8) dx1, dy1, dz1, x(nx+1), y(ny+1), z(nz+1)
    integer i, j, k
    x = 0.d0; y = 0.d0; z = 0.d0  ! Initialize arrays
    dx1 = Lx / dble(nx-6)
    dy1 = Ly / dble(ny-6)
    dz1 = Lz / dble(nz-6)
    dx(:) = dx1
    dy(:) = dy1
    dz(:) = dz1
    ! x direction
    do i = 4, nx-2
      x(i) = dx1 * dble(i-4)
    enddo
    x(1)    = x(4)    - 3.d0 * dx1
    x(2)    = x(4)    - 2.d0 * dx1
    x(3)    = x(4)    - dx1
    x(nx-1) = x(nx-2) + dx1
    x(nx)   = x(nx-2) + 2.d0 * dx1
    x(nx+1) = x(nx-2) + 3.d0 * dx1
    ! y direction
    do j = 4, ny-2
      y(j) = dy1 * dble(j-4)
    enddo
    y(1)    = y(4)    - 3.d0 * dy1
    y(2)    = y(4)    - 2.d0 * dy1
    y(3)    = y(4)    - dy1
    y(ny-1) = y(ny-2) + dy1
    y(ny)   = y(ny-2) + 2.d0 * dy1
    y(ny+1) = y(ny-2) + 3.d0 * dy1
    ! z direction
    do k = 4, nz-2
      z(k) = dz1 * dble(k-4)
    enddo
    z(1)    = z(4)    - 3.d0 * dz1
    z(2)    = z(4)    - 2.d0 * dz1
    z(3)    = z(4)    - dz1
    z(nz-1) = z(nz-2) + dz1
    z(nz)   = z(nz-2) + 2.d0 * dz1
    z(nz+1) = z(nz-2) + 3.d0 * dz1
    ! cell centered
    do i = 1, nx
      xc(i) = 0.5d0 * (x(i) + x(i+1))
    enddo
    do j = 1, ny
      yc(j) = 0.5d0 * (y(j) + y(j+1))
    enddo
    do k = 1, nz
      zc(k) = 0.5d0 * (z(k) + z(k+1))
    enddo
  end subroutine set_grid_cyclic6_3D
end module set_coordinate

