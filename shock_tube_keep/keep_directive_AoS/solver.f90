module solver
  use cudafor
  use parameters, only : dp, Rgas, gamma, Cp, Pr, mu0, T0, S, rhol, Tlr, a, dx, dtdx, nt
  implicit none 
contains 
  !-------------------------------------------------------------------------------------------
  !> set initial condition for Sod shock tube
  subroutine init(nx, q)
    integer, intent(in)     :: nx
    real(dp), intent(inout) :: q(3,nx)
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
      q(1:3,j) = (/rho, rho * u, en/)
    enddo
  end subroutine init

  !-------------------------------------------------------------------------------------------
  !> calculate E (convection term) and Ev (viscous term)
  subroutine eev(nx, q, e, ev)
    integer, intent(in)             :: nx
    real(dp), intent(inout), device :: q(3,nx)
    real(dp), intent(inout), device :: e(3,nx-1), ev(3,nx-1)
    real(dp), device :: qq(3,nx), T(nx), mu(nx)
    real(dp) c, m, g, ie, ke, Lp, tw, twu, kT
    integer ::  j
    ! calculate rho, u, p, T, nu
    !$cuf kernel do(1)<<<32,4>>>
    do j =1, nx
      qq(1,j) = q(1,j)
      qq(2,j) = q(2,j) / qq(1,j)
      qq(3,j) = (gamma - 1.d0) * (q(3,j) - 0.5d0 * qq(1,j) * qq(2,j)**2) ! d0 is unnecessary for u(j)**2
      T(j)    = qq(3,j) / (qq(1,j) * Rgas)
      mu(j)   = mu0 * ((T0 + S) / (T(j) + S)) * (T(j) / T0)**1.5d0
    enddo
    ! calculate KEEP scheme and viscous term
    !$cuf kernel do(1)<<<32,4>>>
    do j = 1, nx-1
      ! KEEP scheme
      ! mass
      c  = 0.25d0 * (qq(1,j) + qq(1,j+1)) * (qq(2,j) + qq(2,j+1))
      ! momentum
      m  = 0.5d0 * c * (qq(2,j) + qq(2,j+1))
      ! pressure gradient
      g  = 0.5d0 * (qq(3,j) + qq(3,j+1))
      ! internal energy
      ie = c * 0.5d0 * Rgas * (T(j) + T(j+1)) / (gamma - 1.d0)
      ! kinetic energy
      ke = 0.5d0 * c * qq(2,j) * qq(2,j+1)
      ! pressure dilatation
      Lp = 0.5d0 * (qq(2,j) * qq(3,j+1) + qq(2,j+1) * qq(3,j))
       
      ! viscous term (2nd-order)
      tw  = 4d0 * (mu(j) + mu(j+1)) * (-qq(2,j) + qq(2,j+1)) / (6d0 * dx)
      twu = 0.5d0 * tw * (qq(2,j) + qq(2,j+1))
      kT  = Cp * (mu(j) + mu(j+1)) * (-T(j) + T(j+1)) / (Pr * 2d0 * dx)

      ! convection term
      e(1,j) = c
      e(2,j) = m + g
      e(3,j) = ie + ke + Lp
      ! viscous term
      ev(1,j) = 0.d0
      ev(2,j) = tw
      ev(3,j) = twu + kT
    enddo
  end subroutine eev

  !-------------------------------------------------------------------------------------------------
  !> 1st-order Euler
  subroutine time_step(nx, q)
    integer, intent(in)             :: nx
    real(dp), intent(inout), device :: q(3,nx)
    real(dp), allocatable, device :: e(:,:), ev(:,:)
    integer :: n, j
    ! cudaMalloc
    allocate(e(3,nx-1), ev(3,nx-1))
    do n = 1, nt
      call  eev(nx, q, e, ev)
      !$cuf kernel do(1)<<<32,4>>>
      do j = 2, nx-1
        q(1:3,j) = q(1:3,j) - dtdx * (-e(1:3,j-1) + e(1:3,j) - (-ev(1:3,j-1) + ev(1:3,j)))
      enddo
      ! boundary conditions are unnecessary in this case
    enddo
    deallocate(e, ev)
  end subroutine time_step
end module solver

