module calc_time_dev
  use mod_globals, only : nt, np
  use calc_flux_base
  use calc_steps
  use set
  use preprocess
  use print
  implicit none
contains
  subroutine RungeKutta(nx, ny, nz, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q)
    integer, intent(in)    :: nx, ny, nz
    real(8), intent(in)    :: x(nx), dx_cpu(nx-1)
    real(8), intent(in)    :: y(ny), dy_cpu(ny-1)
    real(8), intent(in)    :: z(nz), dz_cpu(nz-1), Jacobian_cpu(nx,ny)
    real(8), intent(inout) :: Q(5,nx,ny,nz)
    integer t1, t2, overlap
    ! CPU arrays
    real(8), allocatable :: ruvwp(:,:,:,:), QJ(:,:,:,:), QJs(:,:,:,:), Rs(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable :: T(:,:,:)
    real(8), allocatable :: dx(:), dy(:), dz(:), xix(:), etay(:), zetaz(:), Jacobian(:,:)
    ! for plot
    real(4) :: ke0 = 1.d0, entropy0 = 1.d0

    ! Allocate CPU memory
    allocate(ruvwp(5,nx,ny,nz), T(nx,ny,nz))
    allocate(E(5,nx,ny,nz), F(5,nx,ny,nz), G(5,nx,ny,nz))
    allocate(QJ(5,nx,ny,nz), QJs(5,nx,ny,nz), Rs(5,nx-2,ny-2,nz-2))
    allocate(dx(nx-1), dy(ny-1), dz(nz-1))
    allocate(xix(nx), etay(ny), zetaz(nz), Jacobian(nx,ny))

    print *, "Memory allocation has completed"

    call pre_calc(nx, ny, nz, 0, 1, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q, overlap, &
                  dx, dy, dz, xix, etay, zetaz, Jacobian, QJ, ke0, entropy0)
    Rs = 0.d0

    do t2 = 1, np
      do t1 = 1, nt
        call calc_EFG(nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, ruvwp, T, E, F, G)
        call calc_step(nx, ny, nz, 0.5d0, 1.d0, dx, dy, dz, E, F, G, QJ, QJs, Rs) ! QJs = Q2
        call set_bc(0, nx, ny, nz, Jacobian, QJs)

        call calc_EFG(nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, E, F, G)
        call calc_step(nx, ny, nz, 0.5d0, 2.d0, dx, dy, dz, E, F, G, QJ, QJs, Rs) ! QJs = Q3
        call set_bc(0, nx, ny, nz, Jacobian, QJs)

        call calc_EFG(nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, E, F, G)
        call calc_step(nx, ny, nz, 1.0d0, 2.d0, dx, dy, dz, E, F, G, QJ, QJs, Rs) ! QJs = Q4
        call set_bc(0, nx, ny, nz, Jacobian, QJs)

        call calc_EFG(nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, E, F, G)
        call calc_step4(nx, ny, nz, dx, dy, dz, E, F, G, Rs, QJ)
        call set_bc(0, nx, ny, nz, Jacobian, QJ)
      enddo
      call send_recv_for_print_even(0, 1, t2, nx, ny, nz, x, y, z, Jacobian_cpu, QJ, Q, ke0, entropy0)
    enddo

    deallocate(ruvwp, T, QJ, QJs, Rs, E, F, G, dx, dy, dz, xix, etay, zetaz, Jacobian)
    print *, "Deallocate CPU memory"
  end subroutine RungeKutta
end module calc_time_dev

