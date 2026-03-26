module calc_time_dev
  use mpi
  use mod_globals, only : nt, np
  use calc_flux_base
  use calc_steps
  use set
  use preprocess
  use print
  implicit none
contains
  subroutine RungeKutta(myrank, mygpu, nx, ny, nz, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q)
    integer, intent(in)    :: myrank, mygpu, nx, ny, nz
    real(8), intent(in)    :: x(nx), dx_cpu(nx-1)
    real(8), intent(in)    :: y(ny), dy_cpu(ny-1)
    real(8), intent(in)    :: z(nz), dz_cpu(nz-1), Jacobian_cpu(nx,ny)
    real(8), intent(inout) :: Q(5,nx,ny,nz)
    integer i, j, k, l, t1, t2, overlap, ierr, nranks, stat, ireq, ireq2(2)
    integer istat(MPI_STATUS_SIZE), istat2(MPI_STATUS_SIZE,2)
    real(8), allocatable :: ruvwp(:,:,:,:), QJ(:,:,:,:), QJs(:,:,:,:), Rs(:,:,:,:), E(:,:,:,:), F(:,:,:,:), G(:,:,:,:)
    real(8), allocatable :: T(:,:,:)
    real(8), allocatable :: dx(:), dy(:), dz(:), xix(:), etay(:), zetaz(:), Jacobian(:,:)
    ! for plot
    real(4) :: ke0 = 1.d0, entropy0 = 1.d0

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
    if (myrank == 0) then
      print '(2x, a)', "Running on CPU"
    endif
    if (mod(myrank,2) == 0) then
      call allocate_mem(myrank, nx, ny, nz, dx, dy, dz, xix, etay, zetaz, Jacobian, ruvwp, T, E, F, G)
      allocate(QJ(5,nx,ny,nz), QJs(5,nx,ny,nz), Rs(5,nx-2,ny-2,nz-2))
      print *, "myrank is ", myrank, " memory allocation has completed"
      call pre_calc(nx, ny, nz, myrank, nranks, x, dx_cpu, y, dy_cpu, z, dz_cpu, Jacobian_cpu, Q, overlap, &
                    dx, dy, dz, xix, etay, zetaz, Jacobian, QJ, ke0, entropy0)
      Rs = 0.d0
    else
      call MPI_RECV(ke0,      1, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
      call MPI_RECV(entropy0, 1, MPI_REAL4, myrank-1, myrank,   MPI_COMM_WORLD, istat, ierr)
    endif

    do t2 = 1, np
      do t1 = 1, nt
        if (mod(myrank,2) == 0) then
          call calc_EFG(nx, ny, nz, xix, etay, zetaz, Jacobian, QJ, ruvwp, T, E, F, G)
          call calc_step(nx, ny, nz, 0.5d0, 1.d0, dx, dy, dz, E, F, G, QJ, QJs, Rs)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJs)

          call calc_EFG(nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, E, F, G)
          call calc_step(nx, ny, nz, 0.5d0, 2.d0, dx, dy, dz, E, F, G, QJ, QJs, Rs)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJs)

          call calc_EFG(nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, E, F, G)
          call calc_step(nx, ny, nz, 1.0d0, 2.d0, dx, dy, dz, E, F, G, QJ, QJs, Rs)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJs)

          call calc_EFG(nx, ny, nz, xix, etay, zetaz, Jacobian, QJs, ruvwp, T, E, F, G)
          call calc_step4(nx, ny, nz, dx, dy, dz, E, F, G, Rs, QJ)
          call set_bc(myrank, nx, ny, nz, Jacobian, QJ)
        endif
      enddo
      if (mod(myrank, 2) == 0) then
        call send_recv_for_print_even(myrank, nranks, t2, nx, ny, nz, x, y, z, Jacobian_cpu, QJ, Q, ke0, entropy0)
      else
        call send_recv_for_print_odd(myrank, nranks, t2, nx, ny, nz, x, y, z, Jacobian_cpu, Q, ke0, entropy0)
      endif
    enddo

    if (mod(myrank,2) == 0) then
      deallocate(ruvwp, T, QJ, QJs, Rs, E, F, G, dx, dy, dz, xix, etay, zetaz, Jacobian)
    endif
    print *, "myrank is ", myrank, " deallocate memory"
  end subroutine RungeKutta
end module calc_time_dev
