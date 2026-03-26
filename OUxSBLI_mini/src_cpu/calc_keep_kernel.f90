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
    integer(kind=8), intent(in)    :: id_accuracy
    integer, intent(in)            :: nx, ny, nz
    real(8), intent(in)            :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out)           :: E(5,nx-1,ny-2,nz-2)
    real(4), intent(in), optional  :: sigma(nx,ny,nz)
    integer i, j, k
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j,k)
    do k = 2, nz-1
      do j = 2, ny-1
        do i = 1, nx-1
          if (3 <= i .and. i <= nx-3) then
            E(:,i,j-1,k-1) = KEEP6(Q(1,i-2:i+3,j,k), Q(2,i-2:i+3,j,k), &
                                     Q(3,i-2:i+3,j,k), Q(4,i-2:i+3,j,k), &
                                     Q(2,i-2:i+3,j,k), Q(5,i-2:i+3,j,k), &
                                     T(i-2:i+3,j,k), Normal_x)
          elseif (2 <= i .and. i <= nx-2) then
            E(:,i,j-1,k-1) = KEEP4(Q(1,i-1:i+2,j,k), Q(2,i-1:i+2,j,k), &
                                     Q(3,i-1:i+2,j,k), Q(4,i-1:i+2,j,k), &
                                     Q(2,i-1:i+2,j,k), Q(5,i-1:i+2,j,k), &
                                     T(i-1:i+2,j,k), Normal_x)
          else
            E(:,i,j-1,k-1) = KEEP2(Q(1,i:i+1,j,k), Q(2,i:i+1,j,k), &
                                     Q(3,i:i+1,j,k), Q(4,i:i+1,j,k), &
                                     Q(2,i:i+1,j,k), Q(5,i:i+1,j,k), &
                                     T(i:i+1,j,k), Normal_x)
          endif
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine calc_keep_x6


  subroutine calc_keep_y6(id_accuracy, nx, ny, nz, Q, T, F, sigma)
    use mod_constant, only : Normal_y
    integer(kind=8), intent(in)    :: id_accuracy
    integer, intent(in)            :: nx, ny, nz
    real(8), intent(in)            :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out)           :: F(5,nx-2,ny-1,nz-2)
    real(4), intent(in), optional  :: sigma(nx,ny,nz)
    integer i, j, k
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j,k)
    do k = 2, nz-1
      do j = 1, ny-1
        do i = 2, nx-1
          if (3 <= j .and. j <= ny-3) then
            F(:,i-1,j,k-1) = KEEP6(Q(1,i,j-2:j+3,k), Q(2,i,j-2:j+3,k), &
                                     Q(3,i,j-2:j+3,k), Q(4,i,j-2:j+3,k), &
                                     Q(3,i,j-2:j+3,k), Q(5,i,j-2:j+3,k), &
                                     T(i,j-2:j+3,k), Normal_y)
          elseif (2 <= j .and. j <= ny-2) then
            F(:,i-1,j,k-1) = KEEP4(Q(1,i,j-1:j+2,k), Q(2,i,j-1:j+2,k), &
                                     Q(3,i,j-1:j+2,k), Q(4,i,j-1:j+2,k), &
                                     Q(3,i,j-1:j+2,k), Q(5,i,j-1:j+2,k), &
                                     T(i,j-1:j+2,k), Normal_y)
          else
            F(:,i-1,j,k-1) = KEEP2(Q(1,i,j:j+1,k), Q(2,i,j:j+1,k), &
                                     Q(3,i,j:j+1,k), Q(4,i,j:j+1,k), &
                                     Q(3,i,j:j+1,k), Q(5,i,j:j+1,k), &
                                     T(i,j:j+1,k), Normal_y)
          endif
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine calc_keep_y6


  subroutine calc_keep_z6(id_accuracy, nx, ny, nz, Q, T, G, sigma)
    use mod_constant, only : Normal_z
    integer(kind=8), intent(in)    :: id_accuracy
    integer, intent(in)            :: nx, ny, nz
    real(8), intent(in)            :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out)           :: G(5,nx-2,ny-2,nz-1)
    real(4), intent(in), optional  :: sigma(nx,ny,nz)
    integer i, j, k
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j,k)
    do k = 1, nz-1
      do j = 2, ny-1
        do i = 2, nx-1
          if (3 <= k .and. k <= nz-3) then
            G(:,i-1,j-1,k) = KEEP6(Q(1,i,j,k-2:k+3), Q(2,i,j,k-2:k+3), &
                                     Q(3,i,j,k-2:k+3), Q(4,i,j,k-2:k+3), &
                                     Q(4,i,j,k-2:k+3), Q(5,i,j,k-2:k+3), &
                                     T(i,j,k-2:k+3), Normal_z)
          elseif (2 <= k .and. k <= nz-2) then
            G(:,i-1,j-1,k) = KEEP4(Q(1,i,j,k-1:k+2), Q(2,i,j,k-1:k+2), &
                                     Q(3,i,j,k-1:k+2), Q(4,i,j,k-1:k+2), &
                                     Q(4,i,j,k-1:k+2), Q(5,i,j,k-1:k+2), &
                                     T(i,j,k-1:k+2), Normal_z)
          else
            G(:,i-1,j-1,k) = KEEP2(Q(1,i,j,k:k+1), Q(2,i,j,k:k+1), &
                                     Q(3,i,j,k:k+1), Q(4,i,j,k:k+1), &
                                     Q(4,i,j,k:k+1), Q(5,i,j,k:k+1), &
                                     T(i,j,k:k+1), Normal_z)
          endif
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine calc_keep_z6


  subroutine calc_keep_x4(id_accuracy, nx, ny, nz, Q, T, E, sigma)
    use mod_constant, only : Normal_x
    integer(kind=4), intent(in)    :: id_accuracy
    integer, intent(in)            :: nx, ny, nz
    real(8), intent(in)            :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out)           :: E(5,nx-1,ny-2,nz-2)
    real(4), intent(in), optional  :: sigma(nx,ny,nz)
    integer i, j, k
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j,k)
    do k = 2, nz-1
      do j = 2, ny-1
        do i = 1, nx-1
          if (2 <= i .and. i <= nx-2) then
            E(:,i,j-1,k-1) = KEEP4(Q(1,i-1:i+2,j,k), Q(2,i-1:i+2,j,k), &
                                     Q(3,i-1:i+2,j,k), Q(4,i-1:i+2,j,k), &
                                     Q(2,i-1:i+2,j,k), Q(5,i-1:i+2,j,k), &
                                     T(i-1:i+2,j,k), Normal_x)
          else
            E(:,i,j-1,k-1) = KEEP2(Q(1,i:i+1,j,k), Q(2,i:i+1,j,k), &
                                     Q(3,i:i+1,j,k), Q(4,i:i+1,j,k), &
                                     Q(2,i:i+1,j,k), Q(5,i:i+1,j,k), &
                                     T(i:i+1,j,k), Normal_x)
          endif
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine calc_keep_x4


  subroutine calc_keep_y4(id_accuracy, nx, ny, nz, Q, T, F, sigma)
    use mod_constant, only : Normal_y
    integer(kind=4), intent(in)    :: id_accuracy
    integer, intent(in)            :: nx, ny, nz
    real(8), intent(in)            :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out)           :: F(5,nx-2,ny-1,nz-2)
    real(4), intent(in), optional  :: sigma(nx,ny,nz)
    integer i, j, k
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j,k)
    do k = 2, nz-1
      do j = 1, ny-1
        do i = 2, nx-1
          if (2 <= j .and. j <= ny-2) then
            F(:,i-1,j,k-1) = KEEP4(Q(1,i,j-1:j+2,k), Q(2,i,j-1:j+2,k), &
                                     Q(3,i,j-1:j+2,k), Q(4,i,j-1:j+2,k), &
                                     Q(3,i,j-1:j+2,k), Q(5,i,j-1:j+2,k), &
                                     T(i,j-1:j+2,k), Normal_y)
          else
            F(:,i-1,j,k-1) = KEEP2(Q(1,i,j:j+1,k), Q(2,i,j:j+1,k), &
                                     Q(3,i,j:j+1,k), Q(4,i,j:j+1,k), &
                                     Q(3,i,j:j+1,k), Q(5,i,j:j+1,k), &
                                     T(i,j:j+1,k), Normal_y)
          endif
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine calc_keep_y4


  subroutine calc_keep_z4(id_accuracy, nx, ny, nz, Q, T, G, sigma)
    use mod_constant, only : Normal_z
    integer(kind=4), intent(in)    :: id_accuracy
    integer, intent(in)            :: nx, ny, nz
    real(8), intent(in)            :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out)           :: G(5,nx-2,ny-2,nz-1)
    real(4), intent(in), optional  :: sigma(nx,ny,nz)
    integer i, j, k
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j,k)
    do k = 1, nz-1
      do j = 2, ny-1
        do i = 2, nx-1
          if (2 <= k .and. k <= nz-2) then
            G(:,i-1,j-1,k) = KEEP4(Q(1,i,j,k-1:k+2), Q(2,i,j,k-1:k+2), &
                                     Q(3,i,j,k-1:k+2), Q(4,i,j,k-1:k+2), &
                                     Q(4,i,j,k-1:k+2), Q(5,i,j,k-1:k+2), &
                                     T(i,j,k-1:k+2), Normal_z)
          else
            G(:,i-1,j-1,k) = KEEP2(Q(1,i,j,k:k+1), Q(2,i,j,k:k+1), &
                                     Q(3,i,j,k:k+1), Q(4,i,j,k:k+1), &
                                     Q(4,i,j,k:k+1), Q(5,i,j,k:k+1), &
                                     T(i,j,k:k+1), Normal_z)
          endif
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine calc_keep_z4


  subroutine calc_keep_x2(id_accuracy, nx, ny, nz, Q, T, E)
    use mod_constant, only : Normal_x
    integer(kind=2), intent(in) :: id_accuracy
    integer, intent(in)         :: nx, ny, nz
    real(8), intent(in)         :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out)        :: E(5,nx-1,ny-2,nz-2)
    integer i, j, k
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j,k)
    do k = 2, nz-1
      do j = 2, ny-1
        do i = 1, nx-1
          E(:,i,j-1,k-1) = KEEP2(Q(1,i:i+1,j,k), Q(2,i:i+1,j,k), &
                                   Q(3,i:i+1,j,k), Q(4,i:i+1,j,k), &
                                   Q(2,i:i+1,j,k), Q(5,i:i+1,j,k), &
                                   T(i:i+1,j,k), Normal_x)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine calc_keep_x2


  subroutine calc_keep_y2(id_accuracy, nx, ny, nz, Q, T, F)
    use mod_constant, only : Normal_y
    integer(kind=2), intent(in) :: id_accuracy
    integer, intent(in)         :: nx, ny, nz
    real(8), intent(in)         :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out)        :: F(5,nx-2,ny-1,nz-2)
    integer i, j, k
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j,k)
    do k = 2, nz-1
      do j = 1, ny-1
        do i = 2, nx-1
          F(:,i-1,j,k-1) = KEEP2(Q(1,i,j:j+1,k), Q(2,i,j:j+1,k), &
                                   Q(3,i,j:j+1,k), Q(4,i,j:j+1,k), &
                                   Q(3,i,j:j+1,k), Q(5,i,j:j+1,k), &
                                   T(i,j:j+1,k), Normal_y)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine calc_keep_y2


  subroutine calc_keep_z2(id_accuracy, nx, ny, nz, Q, T, G)
    use mod_constant, only : Normal_z
    integer(kind=2), intent(in) :: id_accuracy
    integer, intent(in)         :: nx, ny, nz
    real(8), intent(in)         :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out)        :: G(5,nx-2,ny-2,nz-1)
    integer i, j, k
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j,k)
    do k = 1, nz-1
      do j = 2, ny-1
        do i = 2, nx-1
          G(:,i-1,j-1,k) = KEEP2(Q(1,i,j,k:k+1), Q(2,i,j,k:k+1), &
                                   Q(3,i,j,k:k+1), Q(4,i,j,k:k+1), &
                                   Q(4,i,j,k:k+1), Q(5,i,j,k:k+1), &
                                   T(i,j,k:k+1), Normal_z)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine calc_keep_z2
end module calc_keep_kernel
