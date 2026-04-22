program main
  use cudafor
  use parameters, only : dp, nx, Rgas, gamma, Lx, rhol, a, Tlr
  use solver
  implicit none
  real(dp) :: q(nx,3)
  real(dp) rho, u, p, x
  real(dp) t_start, t_end
  integer :: n, file_num = 10
  ! gpu memory
  real(dp), allocatable, device :: q_gpu(:,:)

  call init(nx, q)
  call cpu_time(t_start)
  q_gpu = q ! memory transfer cpu -> gpu
  call cpu_time(t_end)
  print *, "Elapsed time (cpu -> gpu) :", t_end - t_start, " [s]"

  call cpu_time(t_start)
  call time_step(nx, q_gpu)
  call cpu_time(t_end)
  print *, "Elapsed time (computation):", t_end - t_start, " [s]"
  
  call cpu_time(t_start)
  q = q_gpu ! memory transfer gpu -> cpu
  call cpu_time(t_end)
  print *, "Elapsed time (gpu -> cpu) :", t_end - t_start, " [s]"
  deallocate(q_gpu)

  open (file_num, file = 'keep.dat')
  do n = 1, nx
    rho = q(n,1)
    u   = q(n,2) / rho
    p   = (gamma - 1.d0) * (q(n,3) - 0.5d0 * rho * u**2)
    x   = dble(n) * dx
    write(file_num,"(4es16.8e3)") x / Lx, rho / rhol, u / a, p / (rhol * Rgas * Tlr)
  enddo 
  close(file_num)
end program main 

