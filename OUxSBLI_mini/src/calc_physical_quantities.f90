module calc_physical_quantities
  use cudafor
  use mod_globals, only : gamma, R
  use mod_constant, only : gamma_1
  implicit none
contains
  subroutine calc_quantities(nx, ny, nz, Jacobian, QJ, Q, T)
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: Jacobian(nx,ny)
    real(8), intent(in), device  :: QJ(5,nx,ny,nz) ! Q / Jacobian
    real(8), intent(out), device :: Q(5,nx,ny,nz), T(nx,ny,nz)
    integer i, j, k
    real(8) :: over_Q1, rho, u, v, w, p
    !$cuf kernel do(3) <<<*,(32,4,2)>>>
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          over_Q1    = 1.d0 / QJ(1,i,j,k)
          rho        = Jacobian(i,j) * QJ(1,i,j,k)
          u          = QJ(2,i,j,k) * over_Q1
          v          = QJ(3,i,j,k) * over_Q1
          w          = QJ(4,i,j,k) * over_Q1
          p          = gamma_1 * (Jacobian(i,j) * QJ(5,i,j,k) - 0.5d0 * rho * (u*u + v*v + w*w))
          Q(1,i,j,k) = rho
          Q(2,i,j,k) = u
          Q(3,i,j,k) = v
          Q(4,i,j,k) = w
          Q(5,i,j,k) = p
          T(i,j,k)   = p / (R * rho)
    enddo;enddo;enddo
  end subroutine calc_quantities
end module calc_physical_quantities

