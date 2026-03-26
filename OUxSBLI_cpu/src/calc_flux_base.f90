module calc_flux_base
  use mod_globals, only : id_accuracy
  use calc_physical_quantities
  use calc_keep_kernel
  implicit none
  private
  public calc_EFG
contains
  subroutine calc_conv_keep(nx, ny, nz, dx, dy, dz, ruvwp, T, E, F, G)
    use mod_globals, only : id_accuracy
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in)  :: dx(nx-1) ! 1 / dx
    real(8), intent(in)  :: dy(ny-1) ! 1 / dy
    real(8), intent(in)  :: dz(nz-1) ! 1 / dz
    real(8), intent(in)  :: ruvwp(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out) :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(out) :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(out) :: G(5,nx-2,ny-2,nz-1)
    call calc_keep_x(id_accuracy, nx, ny, nz, ruvwp, T, E)
    call calc_keep_y(id_accuracy, nx, ny, nz, ruvwp, T, F)
    call calc_keep_z(id_accuracy, nx, ny, nz, ruvwp, T, G)
  end subroutine calc_conv_keep


  subroutine calc_EFG(nx, ny, nz, dx, dy, dz, Jacobian, QJ, ruvwp, T, E, F, G)
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in)  :: dx(nx-1) ! 1 / dx
    real(8), intent(in)  :: dy(ny-1) ! 1 / dy
    real(8), intent(in)  :: dz(nz-1) ! 1 / dz
    real(8), intent(in)  :: Jacobian(nx,ny)
    real(8), intent(in)  :: QJ(5,nx,ny,nz) ! Q / Jacobian
    real(8), intent(out) :: ruvwp(5,nx,ny,nz) ! (rho, u, v, w, p)
    real(8), intent(out) :: T(nx,ny,nz)
    real(8), intent(out) :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(out) :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(out) :: G(5,nx-2,ny-2,nz-1)
    call calc_quantities(nx, ny, nz, Jacobian, QJ, ruvwp, T)
    call calc_conv_keep(nx, ny, nz, dx, dy, dz, ruvwp, T, E, F, G)
  end subroutine calc_EFG
end module calc_flux_base

