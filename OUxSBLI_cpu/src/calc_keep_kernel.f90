module calc_keep_kernel
  use mod_constant, only : R_over_gamma_1, one_third, one_sixth, one_twelfth, two_third
  implicit none
  private
  public calc_keep_x, calc_keep_y, calc_keep_z
  real(8), parameter :: one_24        = 1.d0 / 24.d0
  real(8), parameter :: one_48        = 1.d0 / 48.d0
  real(8), parameter :: one_60        = 1.d0 / 60.d0
  real(8), parameter :: one_120       = 1.d0 / 120.d0
  real(8), parameter :: one_240       = 1.d0 / 240.d0
  real(8), parameter :: seven_twelfth = 7.d0 / 12.d0

  interface calc_keep_x
    module procedure calc_keep_x2, calc_keep_x4, calc_keep_x6
  end interface calc_keep_x

  interface calc_keep_y
    module procedure calc_keep_y2, calc_keep_y4, calc_keep_y6
  end interface calc_keep_y

  interface calc_keep_z
    module procedure calc_keep_z2, calc_keep_z4, calc_keep_z6
  end interface calc_keep_z
contains
  include 'calc_keep_3d.f90'
  subroutine calc_keep_x6(id_accuracy, nx, ny, nz, Q, T, E, sigma)
    use mod_constant, only : Normal_x
    integer(kind=8), intent(in)            :: id_accuracy
    integer, intent(in)                    :: nx, ny, nz
    real(8), intent(in)                    :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out)                   :: E(5,nx-1,ny-2,nz-2)
    real(4), intent(in), optional          :: sigma(nx,ny,nz)
    integer i, j, k
    integer, parameter :: sx = 6 + 5
    real(8), dimension(-1:sx-2) :: rho, u, v, w, p, tmp

    do k = 2, nz-2
      do j = 2, ny-2
        do i = 1, nx-1
          if (i >= 1 .and. i <= nx) then
            rho(-1:sx-2) = Q(1,i-1:i+sx-3,j,k)
            u(-1:sx-2)   = Q(2,i-1:i+sx-3,j,k)
            v(-1:sx-2)   = Q(3,i-1:i+sx-3,j,k)
            w(-1:sx-2)   = Q(4,i-1:i+sx-3,j,k)
            p(-1:sx-2)   = Q(5,i-1:i+sx-3,j,k)
            tmp(-1:sx-2) = T(i-1:i+sx-3,j,k)
          endif

          if (3 <= i .and. i <= nx-3) then
            E(:,i,j-1,k-1) = KEEP6(rho(-1:4), u(-1:4), &
                                     v(-1:4), w(-1:4), &
                                    u(-1:4), p(-1:4), &
                                   tmp(-1:4), Normal_x)
          elseif (2 <= i .and. i <= nx-2) then
            E(:,i,j-1,k-1) = KEEP4(rho(0:3), u(0:3), &
                                     v(0:3), w(0:3), &
                                    u(0:3), p(0:3), &
                                   tmp(0:3), Normal_x)
          else
            E(:,i,j-1,k-1) = KEEP2(rho(0:1), u(0:1), &
                                     v(0:1), w(0:1), &
                                    u(0:1), p(0:1), &
                                   tmp(0:1), Normal_x)
          endif
        enddo
      enddo
    enddo
  end subroutine calc_keep_x6


  subroutine calc_keep_y6(id_accuracy, nx, ny, nz, Q, T, F, sigma)
    use mod_constant, only : Normal_y
    integer(kind=8), intent(in)            :: id_accuracy
    integer, intent(in)                    :: nx, ny, nz
    real(8), intent(in)                    :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out)                   :: F(5,nx-2,ny-1,nz-2)
    real(4), intent(in), optional          :: sigma(nx,ny,nz)
    integer i, j, k
    integer, parameter :: sy = 6 + 5
    real(8), dimension(-1:sy-2) :: rho, u, v, w, p, tmp

    do k = 2, nz-2
      do i = 2, nx-2
        do j = 1, ny-1
          if (j >= 1 .and. j <= ny) then
            rho(-1:sy-2) = Q(1,i,j-1:j+sy-3,k)
            u(-1:sy-2)   = Q(2,i,j-1:j+sy-3,k)
            v(-1:sy-2)   = Q(3,i,j-1:j+sy-3,k)
            w(-1:sy-2)   = Q(4,i,j-1:j+sy-3,k)
            p(-1:sy-2)   = Q(5,i,j-1:j+sy-3,k)
            tmp(-1:sy-2) = T(i,j-1:j+sy-3,k)
          endif

          if (3 <= j .and. j <= ny-3) then
            F(:,i-1,j,k-1) = KEEP6(rho(-1:4), u(-1:4), &
                                     v(-1:4), w(-1:4), &
                                    v(-1:4), p(-1:4), &
                                   tmp(-1:4), Normal_y)
          elseif (2 <= j .and. j <= ny-2) then
            F(:,i-1,j,k-1) = KEEP4(rho(0:3), u(0:3), &
                                     v(0:3), w(0:3), &
                                    v(0:3), p(0:3), &
                                   tmp(0:3), Normal_y)
          else
            F(:,i-1,j,k-1) = KEEP2(rho(0:1), u(0:1), &
                                     v(0:1), w(0:1), &
                                    v(0:1), p(0:1), &
                                   tmp(0:1), Normal_y)
          endif
        enddo
      enddo
    enddo
  end subroutine calc_keep_y6


  subroutine calc_keep_z6(id_accuracy, nx, ny, nz, Q, T, G, sigma)
    use mod_constant, only : Normal_z
    integer(kind=8), intent(in)            :: id_accuracy
    integer, intent(in)                    :: nx, ny, nz
    real(8), intent(in)                    :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out)                   :: G(5,nx-2,ny-2,nz-1)
    real(4), intent(in), optional          :: sigma(nx,ny,nz)
    integer i, j, k
    integer, parameter :: sz = 6 + 5
    real(8), dimension(-1:sz-2) :: rho, u, v, w, p, tmp

    do j = 2, ny-2
      do i = 2, nx-2
        do k = 1, nz-1
          if (k >= 1 .and. k <= nz) then
            rho(-1:sz-2) = Q(1,i,j,k-1:k+sz-3)
            u(-1:sz-2)   = Q(2,i,j,k-1:k+sz-3)
            v(-1:sz-2)   = Q(3,i,j,k-1:k+sz-3)
            w(-1:sz-2)   = Q(4,i,j,k-1:k+sz-3)
            p(-1:sz-2)   = Q(5,i,j,k-1:k+sz-3)
            tmp(-1:sz-2) = T(i,j,k-1:k+sz-3)
          endif

          if (3 <= k .and. k <= nz-3) then
            G(:,i-1,j-1,k) = KEEP6(rho(-1:4), u(-1:4), &
                                     v(-1:4), w(-1:4), &
                                    w(-1:4), p(-1:4), &
                                   tmp(-1:4), Normal_z)
          elseif (2 <= k .and. k <= nz-2) then
            G(:,i-1,j-1,k) = KEEP4(rho(0:3), u(0:3), &
                                     v(0:3), w(0:3), &
                                    w(0:3), p(0:3), &
                                   tmp(0:3), Normal_z)
          else
            G(:,i-1,j-1,k) = KEEP2(rho(0:1), u(0:1), &
                                     v(0:1), w(0:1), &
                                    w(0:1), p(0:1), &
                                   tmp(0:1), Normal_z)
          endif
        enddo
      enddo
    enddo
  end subroutine calc_keep_z6


  subroutine calc_keep_x4(id_accuracy, nx, ny, nz, Q, T, E, sigma)
    use mod_constant, only : Normal_x
    integer(kind=4), intent(in)            :: id_accuracy
    integer, intent(in)                    :: nx, ny, nz
    real(8), intent(in)                    :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out)                   :: E(5,nx-1,ny-2,nz-2)
    real(4), intent(in), optional          :: sigma(nx,ny,nz)
    integer i, j, k
    integer, parameter :: sx = 6 + 3
    real(8), dimension(0:sx-1) :: rho, u, v, w, p, tmp

    do k = 2, nz-2
      do j = 2, ny-2
        do i = 1, nx-1
          if (i >= 1 .and. i <= nx) then
            rho(0:sx-1) = Q(1,i-1:i+sx-2,j,k)
            u(0:sx-1)   = Q(2,i-1:i+sx-2,j,k)
            v(0:sx-1)   = Q(3,i-1:i+sx-2,j,k)
            w(0:sx-1)   = Q(4,i-1:i+sx-2,j,k)
            p(0:sx-1)   = Q(5,i-1:i+sx-2,j,k)
            tmp(0:sx-1) = T(i-1:i+sx-2,j,k)
          endif

          if (2 <= i .and. i <= nx-2) then
            E(:,i,j-1,k-1) = KEEP4(rho(0:3), u(0:3), &
                                     v(0:3), w(0:3), &
                                    u(0:3), p(0:3), &
                                   tmp(0:3), Normal_x)
          else
            E(:,i,j-1,k-1) = KEEP2(rho(0:1), u(0:1), &
                                     v(0:1), w(0:1), &
                                    u(0:1), p(0:1), &
                                   tmp(0:1), Normal_x)
          endif
        enddo
      enddo
    enddo
  end subroutine calc_keep_x4


  subroutine calc_keep_y4(id_accuracy, nx, ny, nz, Q, T, F, sigma)
    use mod_constant, only : Normal_y
    integer(kind=4), intent(in)            :: id_accuracy
    integer, intent(in)                    :: nx, ny, nz
    real(8), intent(in)                    :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out)                   :: F(5,nx-2,ny-1,nz-2)
    real(4), intent(in), optional          :: sigma(nx,ny,nz)
    integer i, j, k
    integer, parameter :: sy = 6 + 3
    real(8), dimension(0:sy-1) :: rho, u, v, w, p, tmp

    do k = 2, nz-2
      do i = 2, nx-2
        do j = 1, ny-1
          if (j >= 1 .and. j <= ny) then
            rho(0:sy-1) = Q(1,i,j-1:j+sy-2,k)
            u(0:sy-1)   = Q(2,i,j-1:j+sy-2,k)
            v(0:sy-1)   = Q(3,i,j-1:j+sy-2,k)
            w(0:sy-1)   = Q(4,i,j-1:j+sy-2,k)
            p(0:sy-1)   = Q(5,i,j-1:j+sy-2,k)
            tmp(0:sy-1) = T(i,j-1:j+sy-2,k)
          endif

          if (2 <= j .and. j <= ny-2) then
            F(:,i-1,j,k-1) = KEEP4(rho(0:3), u(0:3), &
                                     v(0:3), w(0:3), &
                                    v(0:3), p(0:3), &
                                   tmp(0:3), Normal_y)
          else
            F(:,i-1,j,k-1) = KEEP2(rho(0:1), u(0:1), &
                                     v(0:1), w(0:1), &
                                    v(0:1), p(0:1), &
                                   tmp(0:1), Normal_y)
          endif
        enddo
      enddo
    enddo
  end subroutine calc_keep_y4


  subroutine calc_keep_z4(id_accuracy, nx, ny, nz, Q, T, G, sigma)
    use mod_constant, only : Normal_z
    integer(kind=4), intent(in)            :: id_accuracy
    integer, intent(in)                    :: nx, ny, nz
    real(8), intent(in)                    :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out)                   :: G(5,nx-2,ny-2,nz-1)
    real(4), intent(in), optional          :: sigma(nx,ny,nz)
    integer i, j, k
    integer, parameter :: sz = 6 + 3
    real(8), dimension(0:sz-1) :: rho, u, v, w, p, tmp

    do j = 2, ny-2
      do i = 2, nx-2
        do k = 1, nz-1
          if (k >= 1 .and. k <= nz) then
            rho(0:sz-1) = Q(1,i,j,k-1:k+sz-2)
            u(0:sz-1)   = Q(2,i,j,k-1:k+sz-2)
            v(0:sz-1)   = Q(3,i,j,k-1:k+sz-2)
            w(0:sz-1)   = Q(4,i,j,k-1:k+sz-2)
            p(0:sz-1)   = Q(5,i,j,k-1:k+sz-2)
            tmp(0:sz-1) = T(i,j,k-1:k+sz-2)
          endif

          if (2 <= k .and. k <= nz-2) then
            G(:,i-1,j-1,k) = KEEP4(rho(0:3), u(0:3), &
                                     v(0:3), w(0:3), &
                                    w(0:3), p(0:3), &
                                   tmp(0:3), Normal_z)
          else
            G(:,i-1,j-1,k) = KEEP2(rho(0:1), u(0:1), &
                                     v(0:1), w(0:1), &
                                    w(0:1), p(0:1), &
                                   tmp(0:1), Normal_z)
          endif
        enddo
      enddo
    enddo
  end subroutine calc_keep_z4


  subroutine calc_keep_x2(id_accuracy, nx, ny, nz, Q, T, E)
    use mod_constant, only : Normal_x
    integer(kind=2), intent(in) :: id_accuracy
    integer, intent(in)         :: nx, ny, nz
    real(8), intent(in)         :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out)        :: E(5,nx-1,ny-2,nz-2)
    integer i, j, k
    real(8), dimension(2) :: rho, u, v, w, p, tmp

    do k = 2, nz-2
      do j = 2, ny-2
        do i = 1, nx-1
          rho(1) = Q(1,i,j,k)
          rho(2) = Q(1,i+1,j,k)
          u(1)   = Q(2,i,j,k)
          u(2)   = Q(2,i+1,j,k)
          v(1)   = Q(3,i,j,k)
          v(2)   = Q(3,i+1,j,k)
          w(1)   = Q(4,i,j,k)
          w(2)   = Q(4,i+1,j,k)
          p(1)   = Q(5,i,j,k)
          p(2)   = Q(5,i+1,j,k)
          tmp(1) = T(i,j,k)
          tmp(2) = T(i+1,j,k)
          E(:,i,j-1,k-1) = KEEP2(rho, u, v, w, u, p, tmp, Normal_x)
        enddo
      enddo
    enddo
  end subroutine calc_keep_x2


  subroutine calc_keep_y2(id_accuracy, nx, ny, nz, Q, T, F)
    use mod_constant, only : Normal_y
    integer(kind=2), intent(in) :: id_accuracy
    integer, intent(in)         :: nx, ny, nz
    real(8), intent(in)         :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out)        :: F(5,nx-2,ny-1,nz-2)
    integer i, j, k
    real(8), dimension(2) :: rho, u, v, w, p, tmp

    do k = 2, nz-2
      do i = 2, nx-2
        do j = 1, ny-1
          rho(1) = Q(1,i,j,k)
          rho(2) = Q(1,i,j+1,k)
          u(1)   = Q(2,i,j,k)
          u(2)   = Q(2,i,j+1,k)
          v(1)   = Q(3,i,j,k)
          v(2)   = Q(3,i,j+1,k)
          w(1)   = Q(4,i,j,k)
          w(2)   = Q(4,i,j+1,k)
          p(1)   = Q(5,i,j,k)
          p(2)   = Q(5,i,j+1,k)
          tmp(1) = T(i,j,k)
          tmp(2) = T(i,j+1,k)
          F(:,i-1,j,k-1) = KEEP2(rho, u, v, w, v, p, tmp, Normal_y)
        enddo
      enddo
    enddo
  end subroutine calc_keep_y2


  subroutine calc_keep_z2(id_accuracy, nx, ny, nz, Q, T, G)
    use mod_constant, only : Normal_z
    integer(kind=2), intent(in) :: id_accuracy
    integer, intent(in)         :: nx, ny, nz
    real(8), intent(in)         :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out)        :: G(5,nx-2,ny-2,nz-1)
    integer i, j, k
    real(8), dimension(2) :: rho, u, v, w, p, tmp

    do j = 2, ny-2
      do i = 2, nx-2
        do k = 1, nz-1
          rho(1) = Q(1,i,j,k)
          rho(2) = Q(1,i,j,k+1)
          u(1)   = Q(2,i,j,k)
          u(2)   = Q(2,i,j,k+1)
          v(1)   = Q(3,i,j,k)
          v(2)   = Q(3,i,j,k+1)
          w(1)   = Q(4,i,j,k)
          w(2)   = Q(4,i,j,k+1)
          p(1)   = Q(5,i,j,k)
          p(2)   = Q(5,i,j,k+1)
          tmp(1) = T(i,j,k)
          tmp(2) = T(i,j,k+1)
          G(:,i-1,j-1,k) = KEEP2(rho, u, v, w, w, p, tmp, Normal_z)
        enddo
      enddo
    enddo
  end subroutine calc_keep_z2
end module calc_keep_kernel

