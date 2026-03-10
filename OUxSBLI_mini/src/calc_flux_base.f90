module calc_flux_base
  use mod_globals, only : id_accuracy, &
  & blocks, threads, blocksE, blocksF, blocksG, threadsE, threadsF, threadsG
  use calc_physical_quantities
  use calc_keep_kernel
  use set
  implicit none
  private
  public calc_EFG
contains
  subroutine calc_conv_keep(nx, ny, nz, dx, dy, dz, ruvwp, T, E, F, G)
    use mod_globals, only : id_accuracy
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device  :: ruvwp(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out), device :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(out), device :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(out), device :: G(5,nx-2,ny-2,nz-1)
    call calc_keep_x<<<blocksE,threadsE,1>>>(id_accuracy, nx, ny, nz, ruvwp, T, E)
    call calc_keep_y<<<blocksF,threadsF,2>>>(id_accuracy, nx, ny, nz, ruvwp, T, F)
    call calc_keep_z<<<blocksG,threadsG,3>>>(id_accuracy, nx, ny, nz, ruvwp, T, G)
  end subroutine calc_conv_keep


  subroutine calc_EFG(nx, ny, nz, dx, dy, dz, Jacobian, QJ, ruvwp, T, E, F, G)
    integer, intent(in), value   :: nx, ny, nz
    real(8), intent(in), device  :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device  :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device  :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device  :: Jacobian(nx,ny)
    real(8), intent(in), device  :: QJ(5,nx,ny,nz) ! Q / Jacobian
    real(8), intent(out), device :: ruvwp(5,nx,ny,nz) ! (rho, u, v, w, p)
    real(8), intent(out), device :: T(nx,ny,nz)
    real(8), intent(out), device :: E(5,nx-1,ny-2,nz-2)
    real(8), intent(out), device :: F(5,nx-2,ny-1,nz-2)
    real(8), intent(out), device :: G(5,nx-2,ny-2,nz-1)
    integer stat, i, j, k
    call calc_quantities_3D(nx, ny, nz, Jacobian, QJ, ruvwp, T)
    call calc_conv_keep(nx, ny, nz, dx, dy, dz, ruvwp, T, E, F, G)
    stat = cudaDeviceSynchronize()
  end subroutine calc_EFG
end module calc_flux_base

