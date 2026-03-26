program main
  use, intrinsic :: iso_fortran_env
  use mod_globals, only : id_recal, nx, ny, nz, Lx, Ly, Lz
  use set
  use set_coordinate
  use calc_time_dev
  implicit none
  integer i, j, l, m, s, ios
  real(8) t_start, t_end
  real(8), allocatable :: x(:), dx(:), y(:), dy(:), z(:), dz(:), Jacobian(:,:), Q(:,:,:,:)
  character(len=8) header
  character(len=40) filename
  logical is_sequential

  allocate(Q(5,nx,ny,nz), x(nx), dx(nx-1), y(ny), dy(ny-1), z(nz), dz(nz-1), Jacobian(nx,ny))
  call set_grid(0, nx, ny, nz, Lx, Ly, Lz, x, y, z, dx, dy, dz)
  call set_Jacobian_xy(nx, ny, nz, dx, dy, dz, Jacobian)

  if (kind(id_recal) == 4) then
    write(filename, "(a, i5.5, a)") "recal/Q", 1, ".dat"
    open(10, file=filename, action="read", form="unformatted", access="sequential", status="old", iostat=ios)
    if (ios /= 0) then
      print *, "Error opening file."
      stop
    endif
    read(10, iostat=ios) header
    close(10)
    is_sequential = (ios == 0 .and. header == 'SEQFMT01')
    if (is_sequential) then
      open(10, file=filename, action="read", form="unformatted", access="sequential", status="old")
      read(10) header
      read(10) Q
      print *, "simulation has been restarted. access is sequential"
    else
      open(10, file=filename, action="read", form="unformatted", access="stream", status="old")
      rewind(10)
      read(10) Q
      print *, "simulation has been restarted. access is stream"
    endif
    close(10)
  elseif (kind(id_recal) == 2) then
    write(*,*) "set initial condition"
    call set_init(0, nx, ny, nz, x, y, z, Q)
  else
    write(*,*) "wrong paramater was found"
  endif

  call cpu_time(t_start)
  call RungeKutta(nx, ny, nz, x, dx, y, dy, z, dz, Jacobian, Q)
  call cpu_time(t_end)

  ! calculation time
  if (t_end - t_start <= 60.d0) then
    s = int(t_end - t_start)
    print *, "calculation time:", s, " [sec]"
  else
    m = int(t_end - t_start) / 60
    s = int(t_end - t_start) - 60 * m
    print *, "calculation time:", m, " [min] ", s, " [sec]"
  endif
  ! save data
  do l = 1, nz
    do j = 1, ny
      do i = 1, nx
        do m = 1, 5
          Q(m,i,j,l) = Jacobian(i,j) * Q(m,i,j,l)
  enddo;enddo;enddo;enddo
  call cpu_time(t_start)
  write(filename, "(a, i5.5, a)") "recal/Q", 1, ".dat"
  open(10,file=filename,status="replace",action="write",form="unformatted",access="stream")
  write(10) Q
  close(10)
  call cpu_time(t_end)
  s = int(t_end - t_start)
  print *, "output time:", s, " [sec]"

  deallocate(Q, x, dx, y, dy, z, dz, Jacobian)
end program main

