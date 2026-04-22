module solver
  use cudafor
  use parameters, only : dp, Rgas, gamma, Cp, Pr, mu0, T0, S, rhol, Tlr, a, dx, dtdx, nt
  implicit none 
contains 
  !-------------------------------------------------------------------------------------------
  !> set initial condition for Sod shock tube
  subroutine init(nx, q)
    integer, intent(in)     :: nx
    real(dp), intent(inout) :: q(nx,3)
    real(dp)  rho, u, p, en
    integer :: j
    do j = 1, nx
      if (j <= nx/2)then
        rho = rhol
        u   = 0d0
        p   = rhol * Rgas * Tlr
      else
        rho = 0.125d0 * rhol
        u   = 0d0
        p   = 0.1d0 * rhol * Rgas *Tlr
      endif
      en = p / (gamma - 1.d0)
      q(j,1:3) = (/rho, rho * u, en/)
    enddo
  end subroutine init

  !-------------------------------------------------------------------------------------------
  !> calculate E (convection term) and Ev (viscous term)
  subroutine eev(nx, q, e, ev)
    integer, intent(in)             :: nx
    real(dp), intent(inout), device :: q(nx,3)
    real(dp), intent(inout), device :: e(nx-1,3), ev(nx-1,3)
    real(dp), device :: rho(nx), u(nx), p(nx), T(nx), mu(nx)
    real(dp) c, m, g, ie, ke, Lp, tw, twu, kT
    integer ::  j
    ! calculate rho, u, p, T, nu
    !$cuf kernel do(1)<<<32,4>>>
    do j =1, nx
      rho(j) = q(j,1)
      u(j)   = q(j,2) / rho(j)
      p(j)   = (gamma - 1.d0) * (q(j,3) - 0.5d0 * rho(j) * u(j)**2) ! d0 is unnecessary for u(j)**2
      T(j)   = p(j) / (rho(j) * Rgas)
      mu(j)  = mu0 * ((T0 + S) / (T(j) + S)) * (T(j) / T0)**1.5d0
    enddo
    ! calculate KEEP scheme and viscous term
    !$cuf kernel do(1)<<<32,4>>>
    do j = 1, nx-1
      ! KEEP scheme
      ! mass
      c  = 0.25d0 * (rho(j) + rho(j+1)) * (u(j) + u(j+1))
      ! momentum
      m  = 0.5d0 * c * (u(j) + u(j+1))
      ! pressure gradient
      g  = 0.5d0 * (p(j) + p(j+1))
      ! internal energy
      ie = c * 0.5d0 * Rgas * (T(j) + T(j+1)) / (gamma - 1.d0)
      ! kinetic energy
      ke = 0.5d0 * c * u(j) * u(j+1)
      ! pressure dilatation
      Lp = 0.5d0 * (u(j) * p(j+1) + u(j+1) * p(j))
       
      ! viscous term (2nd-order)
      tw  = 4d0 * (mu(j) + mu(j+1)) * (-u(j) + u(j+1)) / (6d0 * dx)
      twu = 0.5d0 * tw * (u(j) + u(j+1))
      kT  = Cp * (mu(j) + mu(j+1)) * (-T(j) + T(j+1)) / (Pr * 2d0 * dx)

      ! convection term
      e(j,1) = c
      e(j,2) = m + g
      e(j,3) = ie + ke + Lp
      ! viscous term
      ev(j,1) = 0.d0
      ev(j,2) = tw
      ev(j,3) = twu + kT
    enddo
  end subroutine eev

  !-------------------------------------------------------------------------------------------------
  !> 1st-order Euler
  subroutine time_step(nx, q)
    integer, intent(in)             :: nx
    real(dp), intent(inout), device :: q(nx,3)
    real(dp), allocatable, device :: e(:,:), ev(:,:)
    integer :: n, j
    ! cudaMalloc
    allocate(e(nx-1,3), ev(nx-1,3))
    do n = 1, nt
      call  eev(nx, q, e, ev)
      !$cuf kernel do(1)<<<32,4>>>
      do j = 2, nx-1
        q(j, 1:3) = q(j, 1:3) - dtdx * (-e(j-1, 1:3) + e(j, 1:3) - (-ev(j-1, 1:3) + ev(j, 1:3)))
      enddo
      ! boundary conditions are unnecessary in this case
    enddo
    deallocate(e, ev)
  end subroutine time_step
end module solver

