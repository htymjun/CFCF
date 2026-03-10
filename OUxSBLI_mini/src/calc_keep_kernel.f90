module calc_keep_kernel
  use mod_globals, only : id_igr, threadsE, threadsF, threadsG
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
 
  attributes(global) subroutine calc_keep_x6(id_accuracy, nx, ny, nz, Q, T, E, sigma)
    use mod_constant, only : Normal_x
    integer(kind=8), intent(in), value    :: id_accuracy
    integer, intent(in), value            :: nx, ny, nz
    real(8), intent(in), device           :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out), device          :: E(5,nx-1,ny-2,nz-2)
    real(4), intent(in), device, optional :: sigma(nx,ny,nz)
    integer i, j, k, it, jt, kt, ii, i_base, idx, offset_yz
    integer, parameter :: sx = threadsE%x + 5
    integer, parameter :: sy = threadsE%y
    integer, parameter :: sz = threadsE%z
    real(8), dimension(-1:sx*sy*sz-2), shared :: rho, u, v, w, p, tmp
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    i_base = (blockIdx%x-1)*blockDim%x
    offset_yz = (jt-1) * sx + (kt-1) * sx * sy
    do ii = it-2, threadsE%x+3, blockDim%x
      i = i_base + ii
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        idx = ii + offset_yz
        rho(idx) = Q(1,i,j,k)
          u(idx) = Q(2,i,j,k)
          v(idx) = Q(3,i,j,k)
          w(idx) = Q(4,i,j,k)
          p(idx) = Q(5,i,j,k)
        tmp(idx) =   T(i,j,k)
      endif
    enddo
    call syncthreads()
    i = (blockIdx%x-1)*blockDim%x + it
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    idx = it + offset_yz
    associate(uu => u)
    if (3 <= i .and. i <= nx-3) then
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(6)
        sig = dble(sigma(i-2:i+3,j,k))
        E(:,i,j-1,k-1) = KEEP_IGR6(rho(idx-2:idx+3), u(idx-2:idx+3), &
                                     v(idx-2:idx+3), w(idx-2:idx+3), &
                                    uu(idx-2:idx+3), p(idx-2:idx+3), &
                                   tmp(idx-2:idx+3), sig, Normal_x)
        end block
      else
        E(:,i,j-1,k-1) = KEEP6(rho(idx-2:idx+3), u(idx-2:idx+3), &
                                 v(idx-2:idx+3), w(idx-2:idx+3), &
                                uu(idx-2:idx+3), p(idx-2:idx+3), &
                               tmp(idx-2:idx+3), Normal_x)
      endif
    elseif (2 <= i .and. i <= nx-2) then
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(4)
        sig = dble(sigma(i-1:i+2,j,k))
        E(:,i,j-1,k-1) = KEEP_IGR4(rho(idx-1:idx+2), u(idx-1:idx+2), &
                                     v(idx-1:idx+2), w(idx-1:idx+2), &
                                    uu(idx-1:idx+2), p(idx-1:idx+2), &
                                   tmp(idx-1:idx+2), sig, Normal_x)
        end block
      else
        E(:,i,j-1,k-1) = KEEP4(rho(idx-1:idx+2), u(idx-1:idx+2), &
                                 v(idx-1:idx+2), w(idx-1:idx+2), &
                                uu(idx-1:idx+2), p(idx-1:idx+2), &
                               tmp(idx-1:idx+2), Normal_x)
      endif
    else
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(2)
        sig = dble(sigma(i:i+1,j,k))
        E(:,i,j-1,k-1) = KEEP_IGR2(rho(idx:idx+1), u(idx:idx+1), &
                                     v(idx:idx+1), w(idx:idx+1), &
                                    uu(idx:idx+1), p(idx:idx+1), &
                                   tmp(idx:idx+1), sig, Normal_x)
        end block
      else
        E(:,i,j-1,k-1) = KEEP2(rho(idx:idx+1), u(idx:idx+1), &
                                 v(idx:idx+1), w(idx:idx+1), &
                                uu(idx:idx+1), p(idx:idx+1), &
                               tmp(idx:idx+1), Normal_x)
      endif
    endif
    end associate
  end subroutine calc_keep_x6


  attributes(global) subroutine calc_keep_y6(id_accuracy, nx, ny, nz, Q, T, F, sigma)
    use mod_constant, only : Normal_y
    integer(kind=8), intent(in), value    :: id_accuracy
    integer, intent(in), value            :: nx, ny, nz
    real(8), intent(in), device           :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out), device          :: F(5,nx-2,ny-1,nz-2)
    real(4), intent(in), device, optional :: sigma(nx,ny,nz)
    integer i, j, k, it, jt, kt, jj, j_base, idx, offset_xz
    integer, parameter :: sx = threadsF%x
    integer, parameter :: sy = threadsF%y + 5
    integer, parameter :: sz = threadsF%z
    real(8), dimension(-1:sx*sy*sz-2), shared :: rho, u, v, w, p, tmp
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    j_base = (blockIdx%y-1)*blockDim%y
    offset_xz = (it-1) * sy + (kt-1) * sy * sx
    do jj = jt-2, threadsF%y+3, blockDim%y
      j = j_base + jj
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        idx = jj + offset_xz
        rho(idx) = Q(1,i,j,k)
          u(idx) = Q(2,i,j,k)
          v(idx) = Q(3,i,j,k)
          w(idx) = Q(4,i,j,k)
          p(idx) = Q(5,i,j,k)
        tmp(idx) =   T(i,j,k)
      endif
    enddo
    call syncthreads()
    j = (blockIdx%y-1)*blockDim%y + jt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    idx = jt + offset_xz
    associate(vv => v)
    if (3 <= j .and. j <= ny-3) then
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(6)
        sig = dble(sigma(i,j-2:j+3,k))
        F(:,i-1,j,k-1) = KEEP_IGR6(rho(idx-2:idx+3), u(idx-2:idx+3), &
                                     v(idx-2:idx+3), w(idx-2:idx+3), &
                                    vv(idx-2:idx+3), p(idx-2:idx+3), &
                                   tmp(idx-2:idx+3), sig, Normal_y)
        end block
      else
        F(:,i-1,j,k-1) = KEEP6(rho(idx-2:idx+3), u(idx-2:idx+3), &
                                 v(idx-2:idx+3), w(idx-2:idx+3), &
                                vv(idx-2:idx+3), p(idx-2:idx+3), &
                               tmp(idx-2:idx+3), Normal_y)
      endif
    elseif (2 <= j .and. j <= ny-2) then
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(4)
        sig = dble(sigma(i,j-1:j+2,k))
        F(:,i-1,j,k-1) = KEEP_IGR4(rho(idx-1:idx+2), u(idx-1:idx+2), &
                                     v(idx-1:idx+2), w(idx-1:idx+2), &
                                    vv(idx-1:idx+2), p(idx-1:idx+2), &
                                   tmp(idx-1:idx+2), sig, Normal_y)
        end block
      else
        F(:,i-1,j,k-1) = KEEP4(rho(idx-1:idx+2), u(idx-1:idx+2), &
                                 v(idx-1:idx+2), w(idx-1:idx+2), &
                                vv(idx-1:idx+2), p(idx-1:idx+2), &
                               tmp(idx-1:idx+2), Normal_y)
      endif
    else
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(2)
        sig = dble(sigma(i,j:j+1,k))
        F(:,i-1,j,k-1) = KEEP_IGR2(rho(idx:idx+1), u(idx:idx+1), &
                                     v(idx:idx+1), w(idx:idx+1), &
                                    vv(idx:idx+1), p(idx:idx+1), &
                                   tmp(idx:idx+1), sig, Normal_y)
        end block
      else
        F(:,i-1,j,k-1) = KEEP2(rho(idx:idx+1), u(idx:idx+1), &
                                 v(idx:idx+1), w(idx:idx+1), &
                                vv(idx:idx+1), p(idx:idx+1), &
                               tmp(idx:idx+1), Normal_y)
      endif
    endif
    end associate
  end subroutine calc_keep_y6


  attributes(global) subroutine calc_keep_z6(id_accuracy, nx, ny, nz, Q, T, G, sigma)
    use mod_constant, only : Normal_z
    integer(kind=8), intent(in), value    :: id_accuracy
    integer, intent(in), value            :: nx, ny, nz
    real(8), intent(in), device           :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out), device          :: G(5,nx-2,ny-2,nz-1)
    real(4), intent(in), device, optional :: sigma(nx,ny,nz)
    integer i, j, k, it, jt, kt, kk, k_base, idx, offset_xy
    integer, parameter :: sx = threadsG%x
    integer, parameter :: sy = threadsG%y
    integer, parameter :: sz = threadsG%z + 5
    real(8), dimension(-1:sx*sy*sz-2), shared :: rho, u, v, w, p, tmp
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k_base = (blockIdx%z-1)*blockDim%z
    offset_xy = (jt-1) * sz + (it-1) * sz * sy
    do kk = kt-2, threadsG%z+3, blockDim%z
      k = k_base + kk
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        idx = kk + offset_xy
        rho(idx) = Q(1,i,j,k)
          u(idx) = Q(2,i,j,k)
          v(idx) = Q(3,i,j,k)
          w(idx) = Q(4,i,j,k)
          p(idx) = Q(5,i,j,k)
        tmp(idx) =   T(i,j,k)
      endif
    enddo
    call syncthreads()
    k = (blockIdx%z-1)*blockDim%z + kt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    idx = kt + offset_xy
    associate(ww => w)
    if (3 <= k .and. k <= nz-3) then
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(6)
        sig = dble(sigma(i,j,k-2:k+3))
        G(:,i-1,j-1,k) = KEEP_IGR6(rho(idx-2:idx+3), u(idx-2:idx+3), &
                                     v(idx-2:idx+3), w(idx-2:idx+3), &
                                    ww(idx-2:idx+3), p(idx-2:idx+3), &
                                   tmp(idx-2:idx+3), sig, Normal_z)
        end block
      else
        G(:,i-1,j-1,k) = KEEP6(rho(idx-2:idx+3), u(idx-2:idx+3), &
                                 v(idx-2:idx+3), w(idx-2:idx+3), &
                                ww(idx-2:idx+3), p(idx-2:idx+3), &
                               tmp(idx-2:idx+3), Normal_z)
      endif
    elseif (2 <= k .and. k <= nz-2) then
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(4)
        sig = dble(sigma(i,j,k-1:k+2))
        G(:,i-1,j-1,kt) = KEEP_IGR4(rho(idx-1:idx+2), u(idx-1:idx+2), &
                                      v(idx-1:idx+2), w(idx-1:idx+2), &
                                     ww(idx-1:idx+2), p(idx-1:idx+2), &
                                    tmp(idx-1:idx+2), sig, Normal_z)
        end block
      else
        G(:,i-1,j-1,k) = KEEP4(rho(idx-1:idx+2), u(idx-1:idx+2), &
                                 v(idx-1:idx+2), w(idx-1:idx+2), &
                                ww(idx-1:idx+2), p(idx-1:idx+2), &
                               tmp(idx-1:idx+2), Normal_z)
      endif
    else
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(2)
        sig = dble(sigma(i,j,k:k+1))
        G(:,i-1,j-1,k) = KEEP_IGR2(rho(idx:idx+1), u(idx:idx+1), &
                                     v(idx:idx+1), w(idx:idx+1), &
                                    ww(idx:idx+1), p(idx:idx+1), &
                                   tmp(idx:idx+1), sig, Normal_z)
        end block
      else
        G(:,i-1,j-1,k) = KEEP2(rho(idx:idx+1), u(idx:idx+1), &
                                 v(idx:idx+1), w(idx:idx+1), &
                                ww(idx:idx+1), p(idx:idx+1), &
                               tmp(idx:idx+1), Normal_z)
      endif
    endif
    end associate
  end subroutine calc_keep_z6


  attributes(global) subroutine calc_keep_x4(id_accuracy, nx, ny, nz, Q, T, E, sigma)
    use mod_constant, only : Normal_x
    integer(kind=4), intent(in), value    :: id_accuracy
    integer, intent(in), value            :: nx, ny, nz
    real(8), intent(in), device           :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out), device          :: E(5,nx-1,ny-2,nz-2)
    real(4), intent(in), device, optional :: sigma(nx,ny,nz)
    integer i, j, k, it, jt, kt, ii, i_base, idx, offset_yz
    integer, parameter :: sx = threadsE%x + 3
    integer, parameter :: sy = threadsE%y
    integer, parameter :: sz = threadsE%z
    real(8), dimension(0:sx*sy*sz-1), shared :: rho, u, v, w, p, tmp
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    i_base = (blockIdx%x-1)*blockDim%x
    offset_yz = (jt-1) * sx + (kt-1) * sx * sy
    do ii = it-1, threadsE%x+2, blockDim%x
      i = i_base + ii
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        idx = ii + offset_yz
        rho(idx) = Q(1,i,j,k)
          u(idx) = Q(2,i,j,k)
          v(idx) = Q(3,i,j,k)
          w(idx) = Q(4,i,j,k)
          p(idx) = Q(5,i,j,k)
        tmp(idx) =   T(i,j,k)
      endif
    enddo
    call syncthreads()
    i = (blockIdx%x-1)*blockDim%x + it
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    idx = it + offset_yz
    associate(uu => u)
    if (2 <= i .and. i <= nx-2) then
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(4)
        sig = dble(sigma(i-1:i+2,j,k))
        E(:,i,j-1,k-1) = KEEP_IGR4(rho(idx-1:idx+2), u(idx-1:idx+2), &
                                     v(idx-1:idx+2), w(idx-1:idx+2), &
                                    uu(idx-1:idx+2), p(idx-1:idx+2), &
                                   tmp(idx-1:idx+2), sig, Normal_x)
        end block
      else
        E(:,i,j-1,k-1) = KEEP4(rho(idx-1:idx+2), u(idx-1:idx+2), &
                                 v(idx-1:idx+2), w(idx-1:idx+2), &
                                uu(idx-1:idx+2), p(idx-1:idx+2), &
                               tmp(idx-1:idx+2), Normal_x)
      endif
    else
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(2)
        sig = dble(sigma(i:i+1,j,k))
        E(:,i,j-1,k-1) = KEEP_IGR2(rho(idx:idx+1), u(idx:idx+1), &
                                     v(idx:idx+1), w(idx:idx+1), &
                                    uu(idx:idx+1), p(idx:idx+1), &
                                   tmp(idx:idx+1), sig, Normal_x)
        end block
      else
        E(:,i,j-1,k-1) = KEEP2(rho(idx:idx+1), u(idx:idx+1), &
                                 v(idx:idx+1), w(idx:idx+1), &
                                uu(idx:idx+1), p(idx:idx+1), &
                               tmp(idx:idx+1), Normal_x)
      endif
    endif
    end associate
  end subroutine calc_keep_x4


  attributes(global) subroutine calc_keep_y4(id_accuracy, nx, ny, nz, Q, T, F, sigma)
    use mod_constant, only : Normal_y
    integer(kind=4), intent(in), value    :: id_accuracy
    integer, intent(in), value            :: nx, ny, nz
    real(8), intent(in), device           :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out), device          :: F(5,nx-2,ny-1,nz-2)
    real(4), intent(in), device, optional :: sigma(nx,ny,nz)
    integer i, j, k, it, jt, kt, jj, j_base, idx, offset_xz
    integer, parameter :: sx = threadsF%x
    integer, parameter :: sy = threadsF%y + 3
    integer, parameter :: sz = threadsF%z
    real(8), dimension(0:sx*sy*sz-1), shared :: rho, u, v, w, p, tmp
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    j_base = (blockIdx%y-1)*blockDim%y
    offset_xz = (it-1) * sy + (kt-1) * sy * sx
    do jj = jt-1, threadsF%y+2, blockDim%y
      j = j_base + jj
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        idx = jj + offset_xz
        rho(idx) = Q(1,i,j,k)
          u(idx) = Q(2,i,j,k)
          v(idx) = Q(3,i,j,k)
          w(idx) = Q(4,i,j,k)
          p(idx) = Q(5,i,j,k)
        tmp(idx) =   T(i,j,k)
      endif
    enddo
    call syncthreads()
    j = (blockIdx%y-1)*blockDim%y + jt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    idx = jt + offset_xz
    associate(vv => v)
    if (2 <= j .and. j <= ny-2) then
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(4)
        sig = dble(sigma(i,j-1:j+2,k))
        F(:,i-1,j,k-1) = KEEP_IGR4(rho(idx-1:idx+2), u(idx-1:idx+2), &
                                     v(idx-1:idx+2), w(idx-1:idx+2), &
                                    vv(idx-1:idx+2), p(idx-1:idx+2), &
                                   tmp(idx-1:idx+2), sig, Normal_y)
        end block
      else
        F(:,i-1,j,k-1) = KEEP4(rho(idx-1:idx+2), u(idx-1:idx+2), &
                                 v(idx-1:idx+2), w(idx-1:idx+2), &
                                vv(idx-1:idx+2), p(idx-1:idx+2), &
                               tmp(idx-1:idx+2), Normal_y)
      endif
    else
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(2)
        sig = dble(sigma(i,j:j+1,k))
        F(:,i-1,j,k-1) = KEEP_IGR2(rho(idx:idx+1), u(idx:idx+1), &
                                     v(idx:idx+1), w(idx:idx+1), &
                                    vv(idx:idx+1), p(idx:idx+1), &
                                   tmp(idx:idx+1), sig, Normal_y)
        end block
      else
        F(:,i-1,j,k-1) = KEEP2(rho(idx:idx+1), u(idx:idx+1), &
                                 v(idx:idx+1), w(idx:idx+1), &
                                vv(idx:idx+1), p(idx:idx+1), &
                               tmp(idx:idx+1), Normal_y)
      endif
    endif
    end associate
  end subroutine calc_keep_y4


  attributes(global) subroutine calc_keep_z4(id_accuracy, nx, ny, nz, Q, T, G, sigma)
    use mod_constant, only : Normal_z
    integer(kind=4), intent(in), value    :: id_accuracy
    integer, intent(in), value            :: nx, ny, nz
    real(8), intent(in), device           :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out), device          :: G(5,nx-2,ny-2,nz-1)
    real(4), intent(in), device, optional :: sigma(nx,ny,nz)
    integer i, j, k, it, jt, kt, kk, k_base, idx, offset_xy
    integer, parameter :: sx = threadsG%x
    integer, parameter :: sy = threadsG%y
    integer, parameter :: sz = threadsG%z + 3
    real(8), dimension(0:sx*sy*sz-1), shared :: rho, u, v, w, p, tmp
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k_base = (blockIdx%z-1)*blockDim%z
    offset_xy = (jt-1) * sz + (it-1) * sz * sy
    do kk = kt-1, threadsG%z+2, blockDim%z
      k = k_base + kk
      if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) then
        idx = kk + offset_xy
        rho(idx) = Q(1,i,j,k)
          u(idx) = Q(2,i,j,k)
          v(idx) = Q(3,i,j,k)
          w(idx) = Q(4,i,j,k)
          p(idx) = Q(5,i,j,k)
        tmp(idx) =   T(i,j,k)
      endif
    enddo
    call syncthreads()
    k = (blockIdx%z-1)*blockDim%z + kt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    idx = kt + offset_xy
    associate(ww => w)
    if (2 <= k .and. k <= nz-2) then
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(4)
        sig = dble(sigma(i,j,k-1:k+2))
        G(:,i-1,j-1,kt) = KEEP_IGR4(rho(idx-1:idx+2), u(idx-1:idx+2), &
                                      v(idx-1:idx+2), w(idx-1:idx+2), &
                                     ww(idx-1:idx+2), p(idx-1:idx+2), &
                                    tmp(idx-1:idx+2), sig, Normal_z)
        end block
      else
        G(:,i-1,j-1,k) = KEEP4(rho(idx-1:idx+2), u(idx-1:idx+2), &
                                 v(idx-1:idx+2), w(idx-1:idx+2), &
                                ww(idx-1:idx+2), p(idx-1:idx+2), &
                               tmp(idx-1:idx+2), Normal_z)
      endif
    else
      if (kind(id_igr) == 4) then
        block
        real(8), device :: sig(2)
        sig = dble(sigma(i,j,k:k+1))
        G(:,i-1,j-1,k) = KEEP_IGR2(rho(idx:idx+1), u(idx:idx+1), &
                                     v(idx:idx+1), w(idx:idx+1), &
                                    ww(idx:idx+1), p(idx:idx+1), &
                                   tmp(idx:idx+1), sig, Normal_z)
        end block
      else
        G(:,i-1,j-1,k) = KEEP2(rho(idx:idx+1), u(idx:idx+1), &
                                 v(idx:idx+1), w(idx:idx+1), &
                                ww(idx:idx+1), p(idx:idx+1), &
                               tmp(idx:idx+1), Normal_z)
      endif
    endif
    end associate
  end subroutine calc_keep_z4


  attributes(global) subroutine calc_keep_x2(id_accuracy, nx, ny, nz, Q, T, E, sigma)
    use mod_constant, only : Normal_x
    integer(kind=2), intent(in), value    :: id_accuracy
    integer, intent(in), value            :: nx, ny, nz
    real(8), intent(in), device           :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out), device          :: E(5,nx-1,ny-2,nz-2)
    real(4), intent(in), device, optional :: sigma(nx,ny,nz)
    integer i, j, k, it, jt, kt
    real(8), dimension(2), device :: rho, u, v, w, p, tmp
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    rho = Q(1,i:i+1,j,k)
    u   = Q(2,i:i+1,j,k)
    v   = Q(3,i:i+1,j,k)
    w   = Q(4,i:i+1,j,k)
    p   = Q(5,i:i+1,j,k)
    tmp = T(i:i+1,j,k)
    if (kind(id_igr) == 4) then
      block
      real(8), device :: sig(2)
      sig = dble(sigma(i:i+1,j,k))
      E(:,i,j-1,k-1) = KEEP_IGR2(rho, u, v, w, u, p, tmp, sig, Normal_x)
      end block
    else
      E(:,i,j-1,k-1) = KEEP2(rho, u, v, w, u, p, tmp, Normal_x)
    endif
  end subroutine calc_keep_x2


  attributes(global) subroutine calc_keep_y2(id_accuracy, nx, ny, nz, Q, T, F, sigma)
    use mod_constant, only : Normal_y
    integer(kind=2), intent(in), value    :: id_accuracy
    integer, intent(in), value            :: nx, ny, nz
    real(8), intent(in), device           :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out), device          :: F(5,nx-2,ny-1,nz-2)
    real(4), intent(in), device, optional :: sigma(nx,ny,nz)
    integer i, j, k, it, jt, kt
    real(8), dimension(2), device :: rho, u, v, w, p, tmp
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    rho = Q(1,i,j:j+1,k)
    u   = Q(2,i,j:j+1,k)
    v   = Q(3,i,j:j+1,k)
    w   = Q(4,i,j:j+1,k)
    p   = Q(5,i,j:j+1,k)
    tmp = T(i,j:j+1,k)
    if (kind(id_igr) == 4) then
      block
      real(8), device :: sig(2)
      sig = dble(sigma(i,j:j+1,k))
      F(:,i-1,j,k-1) = KEEP_IGR2(rho, u, v, w, v, p, tmp, sig, Normal_y)
      end block
    else
      F(:,i-1,j,k-1) = KEEP2(rho, u, v, w, v, p, tmp, Normal_y)
    endif
  end subroutine calc_keep_y2


  attributes(global) subroutine calc_keep_z2(id_accuracy, nx, ny, nz, Q, T, G, sigma)
    use mod_constant, only : Normal_z
    integer(kind=2), intent(in), value    :: id_accuracy
    integer, intent(in), value            :: nx, ny, nz
    real(8), intent(in), device           :: Q(5,nx,ny,nz), T(nx,ny,nz)
    real(8), intent(out), device          :: G(5,nx-2,ny-2,nz-1)
    real(4), intent(in), device, optional :: sigma(nx,ny,nz)
    integer i, j, k, it, jt, kt
    real(8), dimension(2), device :: rho, u, v, w, p, tmp
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    rho = Q(1,i,j,k:k+1)
    u   = Q(2,i,j,k:k+1)
    v   = Q(3,i,j,k:k+1)
    w   = Q(4,i,j,k:k+1)
    p   = Q(5,i,j,k:k+1)
    tmp = T(i,j,k:k+1)
    if (kind(id_igr) == 4) then
      block
      real(8), device :: sig(2)
      sig = dble(sigma(i,j,k:k+1))
      G(:,i-1,j-1,k) = KEEP_IGR2(rho, u, v, w, w, p, tmp, sig, Normal_z)
      end block
    else
      G(:,i-1,j-1,k) = KEEP2(rho, u, v, w, w, p, tmp, Normal_z)
    endif
  end subroutine calc_keep_z2
end module calc_keep_kernel

