module calc_time_dev
  use mod_globals, only : nt, np
  use calc_flux_base
  use calc_steps
  use set
  use print
  implicit none
contains
  subroutine RungeKutta(nx, ny, nz, x, dx, y, dy, z, dz, Jacobian, QJ)
    integer, intent(in)    :: nx, ny, nz
    real(8), intent(in)    :: x(nx), dx(nx-1), y(ny), dy(ny-1), z(nz), dz(nz-1), Jacobian(nx,ny)
    real(8), intent(inout) :: QJ(5,nx,ny,nz)
    integer i, j, k, l, t1, t2
    ! CPU arrays
    real(8), allocatable :: ruvwp(:,:,:,:), QJs(:,:,:,:), Rs(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable :: T(:,:,:)
    ! for plot
    real(4) :: ke0 = 1.d0, entropy0 = 1.d0

    ! Allocate CPU memory
    allocate(ruvwp(5,nx,ny,nz), T(nx,ny,nz))
    allocate(E(5,nx,ny,nz), F(5,nx,ny,nz), G(5,nx,ny,nz))
    allocate(QJs(5,nx,ny,nz), Rs(5,nx-2,ny-2,nz-2))

    print *, "Memory allocation has completed"

    ! set Q / Jacobian
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          do l = 1, 5
            QJ(l,i,j,k) = QJ(l,i,j,k) / Jacobian(i,j)
    enddo;enddo;enddo;enddo
    call print_vtk_3D(0, nx, ny, nz, x, y, z, Jacobian, QJ, ke0, entropy0)
    Rs = 0.d0

    do t2 = 1, np
      do t1 = 1, nt
        call calc_EFG(nx, ny, nz, Jacobian, QJ, ruvwp, T, E, F, G)
        call calc_step(nx, ny, nz, 0.5d0, 1.d0, dx, dy, dz, E, F, G, QJ, QJs, Rs) ! QJs = Q2
        call set_bc(nx, ny, nz, QJs)

        call calc_EFG(nx, ny, nz, Jacobian, QJs, ruvwp, T, E, F, G)
        call calc_step(nx, ny, nz, 0.5d0, 2.d0, dx, dy, dz, E, F, G, QJ, QJs, Rs) ! QJs = Q3
        call set_bc(nx, ny, nz, QJs)

        call calc_EFG(nx, ny, nz, Jacobian, QJs, ruvwp, T, E, F, G)
        call calc_step(nx, ny, nz, 1.0d0, 2.d0, dx, dy, dz, E, F, G, QJ, QJs, Rs) ! QJs = Q4
        call set_bc(nx, ny, nz, QJs)

        call calc_EFG(nx, ny, nz, Jacobian, QJs, ruvwp, T, E, F, G)
        call calc_step4(nx, ny, nz, dx, dy, dz, E, F, G, Rs, QJ)
        call set_bc(nx, ny, nz, QJ)
      enddo
      call print_vtk_3D(t2, nx, ny, nz, x, y, z, Jacobian, QJ, ke0, entropy0)
    enddo

    deallocate(ruvwp, T, QJs, Rs, E, F, G)
    print *, "Deallocate CPU memory"
  end subroutine RungeKutta
end module calc_time_dev
