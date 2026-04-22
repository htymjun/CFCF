module parameters
  implicit none
  integer, parameter  :: dp    = kind(1.0d0)
  integer, parameter  :: nx    = 4096
  integer, parameter  :: nt    = 2500
  real(dp), parameter :: gamma = 1.4d0                                      
  real(dp), parameter :: Rgas  = 287.03d0
  real(dp), parameter :: Pr    = 0.72d0
  real(dp), parameter :: mu0   = 1.716d-5
  real(dp), parameter :: T0    = 273.2d0
  real(dp), parameter :: S     = 111d0
  real(dp), parameter :: Re    = 25000d0
  real(dp), parameter :: Tlr   = 300.d0
  real(dp), parameter :: CFL   = 0.1d0 
  real(dp), parameter :: rhol  = 1.293d0
  !> don't have to edit
  real(dp), parameter :: Cp    = gamma * Rgas / (gamma - 1.d0)
  real(dp), parameter :: a     = sqrt(Rgas * 300d0)  
  real(dp), parameter :: Lx    = Re * mu0 / (1.293d0 * a)
  real(dp), parameter :: dx    = Lx / dble(nx-1)
  real(dp), parameter :: dtdx  = CFL * dx / (a * dx)
end module parameters

