module preprocess
  use cudafor
  use mpi
  use print
  implicit none
contains
  subroutine check_gpu(mygpu)
    integer, intent(in) :: mygpu
    integer ilen, stat
    type(cudaDeviceProp) prop
    stat = cudaSetDevice(mygpu)
    stat = cudaGetDeviceProperties(prop, mygpu)
    ilen = verify(prop%name, ' ', .true.)
    print '(1x, a, a, i1, a)', prop%name(1:ilen), " (GPU", mygpu, ") is available"
  end subroutine check_gpu


  subroutine allocate_device_mem(myrank, nx, ny, nz, dx, dy, dz, xix, etay, zetaz, Jacobian, ruvwp, T, E, F, G)
    integer, intent(in)                       :: myrank, nx, ny, nz
    real(8), intent(out), allocatable, device :: dx(:), dy(:), dz(:), xix(:), etay(:), zetaz(:), Jacobian(:,:)
    real(8), intent(out), allocatable, device :: ruvwp(:,:,:,:), T(:,:,:)
    real(8), intent(out), allocatable, device :: E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
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


  subroutine pre_calc(nx, ny, nz, myrank, nranks, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q, overlap, &
                      dx, dy, dz, xix, etay, zetaz, Jacobian, QJ, ke0, entropy0)
    integer, intent(in)    :: nx, ny, nz, myrank, nranks
    real(8), intent(in)    :: x(nx), dx_cpu(nx-1), y(ny), dy_cpu(ny-1), z(nz), dz_cpu(nz-1), Jacobian_cpu(nx,ny)
    real(8), intent(inout) :: Q(5,nx,ny,nz)
    integer, intent(out)   :: overlap
    real(8), intent(out), device :: dx(nx-1), dy(ny-1), dz(nz-1), xix(nx-1), etay(ny-1), zetaz(nz-1), Jacobian(nx,ny)
    real(8), intent(out), device :: QJ(5,nx,ny,nz)
    real(4), intent(inout)       :: ke0, entropy0
    real(8) xix_cpu(nx-1), etay_cpu(ny-1), zetaz_cpu(nz-1)
    real(4) rho1d(nx*ny*nz), p1d(nx*ny*nz), v1d(nx*ny*nz*3)
    integer i, j, k, l, ierr
    ! set Q / Jacobian
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          do l = 1, 5
            Q(l,i,j,k) = Q(l,i,j,k) / Jacobian_cpu(i,j)
    enddo;enddo;enddo;enddo
    ! copy on GPU
    xix_cpu   = 1.d0 / dx_cpu
    etay_cpu  = 1.d0 / dy_cpu
    zetaz_cpu = 1.d0 / dz_cpu
    dx       = dx_cpu
    dy       = dy_cpu
    dz       = dz_cpu
    xix      = xix_cpu
    etay     = etay_cpu
    zetaz    = zetaz_cpu
    Jacobian = Jacobian_cpu
    QJ = Q
    ! for multi GPU
    if (kind(id_accuracy) == 8) then
      overlap = 3
    elseif (kind(id_accuracy) == 4) then
      overlap = 2
    else
      overlap = 1
    endif
    call make_1d_for_print(nx, ny, nz, Jacobian_cpu, Q, rho1d, p1d, v1d)
    call print_vtk(0, nx, ny, nz, myrank+1, nranks, x, y, z, rho1d, p1d, v1d, ke0, entropy0)
    call MPI_SEND(ke0,      1, MPI_REAL4, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
    call MPI_SEND(entropy0, 1, MPI_REAL4, myrank+1, myrank+1, MPI_COMM_WORLD, ierr)
  end subroutine pre_calc
end module preprocess

