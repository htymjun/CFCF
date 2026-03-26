module preprocess
  use print
  implicit none
contains
  subroutine check_gpu()
    ! CPU version - no GPU checking needed
  end subroutine check_gpu


  subroutine allocate_device_mem(myrank, nx, ny, nz, dx, dy, dz, xix, etay, zetaz, Jacobian, ruvwp, T, E, F, G)
    integer, intent(in)                       :: myrank, nx, ny, nz
    real(8), intent(out), allocatable :: dx(:), dy(:), dz(:), xix(:), etay(:), zetaz(:), Jacobian(:,:)
    real(8), intent(out), allocatable :: ruvwp(:,:,:,:), T(:,:,:)
    real(8), intent(out), allocatable :: E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    integer ierr
    allocate(ruvwp(5,nx,ny,nz), E(5,nx-1,ny-2,nz-2), F(5,nx-2,ny-1,nz-2), G(5,nx-2,ny-2,nz-1), stat=ierr)
    allocate(dx(nx-1), dy(ny-1), dz(nz-1), xix(nx-1), etay(ny-1), zetaz(nz-1), Jacobian(nx,ny), stat=ierr)
    allocate(T(nx,ny,nz), stat=ierr)
    if (ierr /= 0) then
      print *, "myrank is ", myrank, " memory allocation failed", ierr
    else
      print *, "myrank is ", myrank, " memory allocation has completed"
    endif
  end subroutine allocate_device_mem


  subroutine pre_calc(nx, ny, nz, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q, &
                      dx, dy, dz, Jacobian, QJ, ke0, entropy0)
    integer, intent(in)    :: nx, ny, nz
    real(8), intent(in)    :: x(nx), dx_cpu(nx-1), y(ny), dy_cpu(ny-1), z(nz), dz_cpu(nz-1), Jacobian_cpu(nx,ny)
    real(8), intent(inout) :: Q(5,nx,ny,nz)
    real(8), intent(out)   :: dx(nx-1), dy(ny-1), dz(nz-1), Jacobian(nx,ny)
    real(8), intent(out)   :: QJ(5,nx,ny,nz)
    real(4), intent(inout) :: ke0, entropy0
    integer i, j, k, l
    ! set Q / Jacobian
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          do l = 1, 5
            Q(l,i,j,k) = Q(l,i,j,k) / Jacobian_cpu(i,j)
    enddo;enddo;enddo;enddo
    call print_vtk_3D(0, nx, ny, nz, x, y, z, Jacobian, Q, ke0, entropy0)
  end subroutine pre_calc
end module preprocess

