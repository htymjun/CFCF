module calc_keep_kernel
  use mod_constant, only : R_over_gamma_1
  implicit none
  private
  public calc_keep_x, calc_keep_y, calc_keep_z
contains
  pure function KEEP2(rho, u, v, w, uu, p, T, Normal) result(F)
    real(8), intent(in), dimension(2) :: rho, u, v, w, uu, p, T
    real(8), intent(in), dimension(5) :: Normal
    real(8) F(5)
    F(1) = 0.25d0 * (rho(1) + rho(2)) * (uu(1) + uu(2))
    F(2) = 0.5d0 * (F(1) * (u(1) + u(2)) + (p(1) + p(2)) * Normal(2))
    F(3) = 0.5d0 * (F(1) * (v(1) + v(2)) + (p(1) + p(2)) * Normal(3))
    F(4) = 0.5d0 * (F(1) * (w(1) + w(2)) + (p(1) + p(2)) * Normal(4))
    F(5) = F(1) * 0.5d0 * (T(1) + T(2)) * R_over_gamma_1 ! internal energy
    F(5) = F(5) + 0.5d0 * (uu(1) * p(2) + uu(2) * p(1)) ! pressure diffusion
    F(5) = F(5) + 0.5d0 * F(1) * (u(1) * u(2) + v(1) * v(2) + w(1) * w(2)) ! kinetic energy
  end function KEEP2

  subroutine calc_keep_x(nx, ny, nz, Q, T, E)
    use mod_constant, only : Normal_x
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(in)  :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out) :: E(5,nx-1,ny-2,nz-2)
    integer i, j, k
    real(8), dimension(2) :: rho, u, v, w, p, tmp
    do k = 2, nz-2
      do j = 2, ny-2
        do i = 1, nx-1
          rho(1) = Q(1,i,j,k)
          rho(2) = Q(1,i+1,j,k)
          u(1)   = Q(2,i,j,k)
          u(2)   = Q(2,i+1,j,k)
          v(1)   = Q(3,i,j,k)
          v(2)   = Q(3,i+1,j,k)
          w(1)   = Q(4,i,j,k)
          w(2)   = Q(4,i+1,j,k)
          p(1)   = Q(5,i,j,k)
          p(2)   = Q(5,i+1,j,k)
          tmp(1) = T(i,j,k)
          tmp(2) = T(i+1,j,k)
          E(:,i,j-1,k-1) = KEEP2(rho, u, v, w, u, p, tmp, Normal_x)
    enddo;enddo;enddo
  end subroutine calc_keep_x


  subroutine calc_keep_y(nx, ny, nz, Q, T, F)
    use mod_constant, only : Normal_y
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(in)  :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out) :: F(5,nx-2,ny-1,nz-2)
    integer i, j, k
    real(8), dimension(2) :: rho, u, v, w, p, tmp
    do k = 2, nz-2
      do i = 2, nx-2
        do j = 1, ny-1
          rho(1) = Q(1,i,j,k)
          rho(2) = Q(1,i,j+1,k)
          u(1)   = Q(2,i,j,k)
          u(2)   = Q(2,i,j+1,k)
          v(1)   = Q(3,i,j,k)
          v(2)   = Q(3,i,j+1,k)
          w(1)   = Q(4,i,j,k)
          w(2)   = Q(4,i,j+1,k)
          p(1)   = Q(5,i,j,k)
          p(2)   = Q(5,i,j+1,k)
          tmp(1) = T(i,j,k)
          tmp(2) = T(i,j+1,k)
          F(:,i-1,j,k-1) = KEEP2(rho, u, v, w, v, p, tmp, Normal_y)
    enddo;enddo;enddo
  end subroutine calc_keep_y


  subroutine calc_keep_z(nx, ny, nz, Q, T, G)
    use mod_constant, only : Normal_z
    integer, intent(in)  :: nx, ny, nz
    real(8), intent(in)  :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out) :: G(5,nx-2,ny-2,nz-1)
    integer i, j, k
    real(8), dimension(2) :: rho, u, v, w, p, tmp
    do j = 2, ny-2
      do i = 2, nx-2
        do k = 1, nz-1
          rho(1) = Q(1,i,j,k)
          rho(2) = Q(1,i,j,k+1)
          u(1)   = Q(2,i,j,k)
          u(2)   = Q(2,i,j,k+1)
          v(1)   = Q(3,i,j,k)
          v(2)   = Q(3,i,j,k+1)
          w(1)   = Q(4,i,j,k)
          w(2)   = Q(4,i,j,k+1)
          p(1)   = Q(5,i,j,k)
          p(2)   = Q(5,i,j,k+1)
          tmp(1) = T(i,j,k)
          tmp(2) = T(i,j,k+1)
          G(:,i-1,j-1,k) = KEEP2(rho, u, v, w, w, p, tmp, Normal_z)
    enddo;enddo;enddo
  end subroutine calc_keep_z
end module calc_keep_kernel
