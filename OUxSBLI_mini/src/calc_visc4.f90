module calc_visc4
  use mod_globals, only : id_visc, gamma, R, Pr, Prt, dt, threadsEv, threadsFv, threadsGv
  use mod_constant, only : Cp, gamma_1, Cp_over_Pr, one_third, two_third, one_twelfth
  implicit none
  private
  public calc_Ev4, calc_Fv4, calc_Gv4
  real(8), parameter :: one_24 = 1.d0 / 24.d0
contains
  pure attributes(device) function flux4(a) result(ans)
    real(8), intent(in), device :: a(3)
    real(8) ans
    ans = (-a(1) + 26.d0 * a(2) - a(3)) * one_24
  end function flux4


  pure attributes(device) subroutine calc_tau_straight(mu, u, vy, wz, d, t11, ut11)
    real(8), intent(in), contiguous :: mu(3), u(6), vy(6), wz(6)
    real(8), intent(in)             :: d
    real(8), intent(out)            :: t11, ut11
    real(8) tmp1, tmp2, tmp3
    tmp1 = two_third * mu(1) * ((2.25d0 * (-u(2) + u(3)) - (-u(1) + u(4)) * one_twelfth) * d &
           - 0.0625d0 * (-vy(1) + 9.d0 * (vy(2) + vy(3)) - vy(4)) &
           - 0.0625d0 * (-wz(1) + 9.d0 * (wz(2) + wz(3)) - wz(4)))
    tmp2 = two_third * mu(2) * ((2.25d0 * (-u(3) + u(4)) - (-u(2) + u(5)) * one_twelfth) * d &
           - 0.0625d0 * (-vy(2) + 9.d0 * (vy(3) + vy(4)) - vy(5)) &
           - 0.0625d0 * (-wz(2) + 9.d0 * (wz(3) + wz(4)) - wz(5)))
    tmp3 = two_third * mu(3) * ((2.25d0 * (-u(4) + u(5)) - (-u(3) + u(6)) * one_twelfth) * d &
           - 0.0625d0 * (-vy(3) + 9.d0 * (vy(4) + vy(5)) - vy(6)) &
           - 0.0625d0 * (-wz(3) + 9.d0 * (wz(4) + wz(5)) - wz(6)))
    t11  = (-tmp1 + 26.d0 * tmp2 - tmp3) * one_24
    tmp1 = 0.0625d0 * (-u(1) + 9.d0 * (u(2) + u(3)) - u(4)) * tmp1
    tmp2 = 0.0625d0 * (-u(2) + 9.d0 * (u(3) + u(4)) - u(5)) * tmp2
    tmp3 = 0.0625d0 * (-u(3) + 9.d0 * (u(4) + u(5)) - u(6)) * tmp3
    ut11 = (-tmp1 + 26.d0 * tmp2 - tmp3) * one_24
  end subroutine calc_tau_straight


  pure attributes(device) subroutine calc_tau_cross(mu, v, uy, d, t12, vt12)
    real(8), intent(in), contiguous :: mu(3), v(6), uy(6)
    real(8), intent(in)             :: d
    real(8), intent(out)            :: t12, vt12
    real(8) tmp1, tmp2, tmp3
    tmp1 = mu(1) * ((1.125d0 * (-v(2) + v(3)) - (-v(1) + v(4)) * one_24) * d &
                    + 0.0625d0 * (-uy(1) + 9.d0 * (uy(2) + uy(3)) - uy(4)))
    tmp2 = mu(2) * ((1.125d0 * (-v(3) + v(4)) - (-v(2) + v(5)) * one_24) * d &
                    + 0.0625d0 * (-uy(2) + 9.d0 * (uy(3) + uy(4)) - uy(5)))
    tmp3 = mu(3) * ((1.125d0 * (-v(4) + v(5)) - (-v(3) + v(6)) * one_24) * d &
                    + 0.0625d0 * (-uy(3) + 9.d0 * (uy(4) + uy(5)) - uy(6)))
    t12  = (-tmp1 + 26.d0 * tmp2 - tmp3) * one_24
    tmp1 = 0.0625d0 * (-v(1) + 9.d0 * (v(2) + v(3)) - v(4)) * tmp1
    tmp2 = 0.0625d0 * (-v(2) + 9.d0 * (v(3) + v(4)) - v(5)) * tmp2
    tmp3 = 0.0625d0 * (-v(3) + 9.d0 * (v(4) + v(5)) - v(6)) * tmp3
    vt12 = (-tmp1 + 26.d0 * tmp2 - tmp3) * one_24
  end subroutine calc_tau_cross


  attributes(global) subroutine calc_Ev4(nx, ny, nz, dx, dy, dz, Q, T, mu, E)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz), T(nx,ny,nz), mu(nx,ny,nz)
    real(8), intent(inout), device :: E(5,nx-1,ny-2,nz-2)
    real(8), shared ::  u(-2:threadsEv%x+3,threadsEv%y,threadsEv%z)
    real(8), shared ::  v(-2:threadsEv%x+3,threadsEv%y,threadsEv%z)
    real(8), shared ::  w(-2:threadsEv%x+3,threadsEv%y,threadsEv%z)
    real(8), shared :: uy(-2:threadsEv%x+3,threadsEv%y,threadsEv%z)
    real(8), shared :: vy(-2:threadsEv%x+3,threadsEv%y,threadsEv%z)
    real(8), shared :: uz(-2:threadsEv%x+3,threadsEv%z,threadsEv%y)
    real(8), shared :: wz(-2:threadsEv%x+3,threadsEv%z,threadsEv%y)
    integer i, j, k, it, jt, kt, ii, i_base
    real(8) :: txx, txy, txz, utxx, vtxy, wtxz, kTx
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    i_base = (blockIdx%x-1)*blockDim%x
    do ii = it-2, threadsEv%x+3, blockDim%x
      i = i_base + ii
      if (1 <= i .and. i <= nx .and. j <= ny .and. k <= nz) then
        u(ii,jt,kt) = Q(2,i,j,k)
        v(ii,jt,kt) = Q(3,i,j,k)
        w(ii,jt,kt) = Q(4,i,j,k)
      endif
      if (1 <= i .and. i <= nx .and. 3 <= j .and. j <= ny-2 .and. k <= nz) then
        uy(ii,jt,kt) = (two_third * (-Q(2,i,j-1,k) + Q(2,i,j+1,k)) - one_twelfth * (-Q(2,i,j-2,k) + Q(2,i,j+2,k))) * dy(j)
        vy(ii,jt,kt) = (two_third * (-Q(3,i,j-1,k) + Q(3,i,j+1,k)) - one_twelfth * (-Q(3,i,j-2,k) + Q(3,i,j+2,k))) * dy(j)
      endif
      if (1 <= i .and. i <= nx .and. j <= ny .and. 3 <= k .and. k <= nz-2) then
        uz(ii,jt,kt) = (two_third * (-Q(2,i,j,k-1) + Q(2,i,j,k+1)) - one_twelfth * (-Q(2,i,j,k-2) + Q(2,i,j,k+2))) * dz(k)
        wz(ii,jt,kt) = (two_third * (-Q(4,i,j,k-1) + Q(4,i,j,k+1)) - one_twelfth * (-Q(4,i,j,k-2) + Q(4,i,j,k+2))) * dz(k)
      endif
    enddo
    call syncthreads()
    i  = (blockIdx%x-1)*blockDim%x + it
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    if (3 <= i .and. i <= nx-3 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-2) then
      block
        real(8), device :: mu3(3)
        mu3(:) = 0.0625d0 * (9.d0 * (mu(i-1:i+1,j,k) + mu(i:i+2,j,k)) - (mu(i-2:i,j,k) + mu(i+1:i+3,j,k)))
        block ! dQdx
          real(8), device :: kTx3(3)
          kTx3(:) = Cp_over_Pr * mu3(:) * &
                    (1.125d0 * (-T(i-1:i+1,j,k) + T(i:i+2,j,k)) - (-T(i-2:i,j,k) + T(i+1:i+3,j,k)) * one_24) * dx(i)
          kTx     = flux4(kTx3(:))
        end block
        call calc_tau_straight(mu3, u(it-2:it+3,jt,kt), vy(it-2:it+3,jt,kt), wz(it-2:it+3,jt,kt), dx(i), txx, utxx)
        call calc_tau_cross(mu3, v(it-2:it+3,jt,kt), uy(it-2:it+3,jt,kt), dx(i), txy, vtxy)
        call calc_tau_cross(mu3, w(it-2:it+3,jt,kt), uz(it-2:it+3,jt,kt), dx(i), txz, wtxz)
      end block
    else
      block
        real(8) mx, mux, mvx, mwx, muy, mvy, muz, mwz
        mx  = 0.5d0 * (mu(i,j,k) + mu(i+1,j,k))
        kTx = Cp_over_Pr * mx * (-T(i,j,k) + T(i+1,j,k)) * dx(i)
        block
          real(8), device :: my(2)
          my(1) = 0.25d0 * (mu(i,j-1,k) + mu(i,j,  k) + mu(i+1,j-1,k) + mu(i+1,j,  k))
          my(2) = 0.25d0 * (mu(i,j,  k) + mu(i,j+1,k) + mu(i+1,j,  k) + mu(i+1,j+1,k))
          muy = 0.25d0 * (my(1) * (-Q(2,i,j-1,k) + Q(2,i,j,k) - Q(2,i+1,j-1,k) + Q(2,i+1,j,k)) &
                        + my(2) * (-Q(2,i,j,k) + Q(2,i,j+1,k) - Q(2,i+1,j,k) + Q(2,i+1,j+1,k))) * dy(j)
          mvy = 0.25d0 * (my(1) * (-Q(3,i,j-1,k) + Q(3,i,j,k) - Q(3,i+1,j-1,k) + Q(3,i+1,j,k)) &
                        + my(2) * (-Q(3,i,j,k) + Q(3,i,j+1,k) - Q(3,i+1,j,k) + Q(3,i+1,j+1,k))) * dy(j)
        end block
        block
          real(8), device :: mz(2)
          mz(1) = 0.25d0 * (mu(i,j,k-1) + mu(i,j,k  ) + mu(i+1,j,k-1) + mu(i+1,j,k  ))
          mz(2) = 0.25d0 * (mu(i,j,k  ) + mu(i,j,k+1) + mu(i+1,j,k  ) + mu(i+1,j,k+1))
          muz = 0.25d0 * (mz(1) * (-Q(2,i,j,k-1) + Q(2,i,j,k) - Q(2,i+1,j,k-1) + Q(2,i+1,j,k)) &
                        + mz(2) * (-Q(2,i,j,k) + Q(2,i,j,k+1) - Q(2,i+1,j,k) + Q(2,i+1,j,k+1))) * dz(k)
          mwz = 0.25d0 * (mz(1) * (-Q(4,i,j,k-1) + Q(4,i,j,k) - Q(4,i+1,j,k-1) + Q(4,i+1,j,k)) &
                        + mz(2) * (-Q(4,i,j,k) + Q(4,i,j,k+1) - Q(4,i+1,j,k) + Q(4,i+1,j,k+1))) * dz(k)
        end block
        mux  = mx * (-u(it,jt,kt) + u(it+1,jt,kt)) * dx(i)
        mvx  = mx * (-v(it,jt,kt) + v(it+1,jt,kt)) * dx(i)
        mwx  = mx * (-v(it,jt,kt) + w(it+1,jt,kt)) * dx(i)
        txx  = two_third * (2.d0 * mux - mvy - mwz)
        txy  = muy + mvx
        txz  = mwx + muz
        utxx = 0.5d0 * (u(it,jt,kt) + u(it+1,jt,kt)) * txx
        vtxy = 0.5d0 * (v(it,jt,kt) + v(it+1,jt,kt)) * txy
        wtxz = 0.5d0 * (w(it,jt,kt) + w(it+1,jt,kt)) * txz
      end block
    endif
    E(2,i,j-1,k-1) = E(2,i,j-1,k-1) - txx
    E(3,i,j-1,k-1) = E(3,i,j-1,k-1) - txy
    E(4,i,j-1,k-1) = E(4,i,j-1,k-1) - txz
    E(5,i,j-1,k-1) = E(5,i,j-1,k-1) - (utxx + vtxy + wtxz + kTx)
  end subroutine calc_Ev4
 

  attributes(global) subroutine calc_Fv4(nx, ny, nz, dy, dx, dz, Q, T, mu, F)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz), T(nx,ny,nz), mu(nx,ny,nz)
    real(8), intent(inout), device :: F(5,nx-2,ny-1,nz-2)
    real(8), shared ::  u(-2:threadsFv%y+3,threadsFv%x,threadsFv%z)
    real(8), shared ::  v(-2:threadsFv%y+3,threadsFv%x,threadsFv%z)
    real(8), shared ::  w(-2:threadsFv%y+3,threadsFv%x,threadsFv%z)
    real(8), shared :: ux(-2:threadsFv%y+3,threadsFv%x,threadsFv%z)
    real(8), shared :: vx(-2:threadsFv%y+3,threadsFv%x,threadsFv%z)
    real(8), shared :: vz(-2:threadsFv%y+3,threadsFv%x,threadsFv%z)
    real(8), shared :: wz(-2:threadsFv%y+3,threadsFv%x,threadsFv%z)
    integer i, j, k, it, jt, kt, jj, j_base
    real(8) :: tyx, tyy, tyz, utyx, vtyy, wtyz, kTy
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    j_base = (blockIdx%y-1)*blockDim%y
    do jj = jt-2, threadsFv%y+3, blockDim%y
      j = j_base + jj
      if (i <= nx .and. 1 <= j .and. j <= ny .and. k <= nz) then
        u(jj,it,kt) = Q(2,i,j,k)
        v(jj,it,kt) = Q(3,i,j,k)
        w(jj,it,kt) = Q(4,i,j,k)
      endif
      if (3 <= i .and. i <= nx-2 .and. 1 <= j .and. j <= ny .and. k <= nz) then
        ux(jj,it,kt) = (two_third * (-Q(2,i-1,j,k) + Q(2,i+1,j,k)) - one_twelfth * (-Q(2,i-2,j,k) + Q(2,i+2,j,k))) * dx(i)
        vx(jj,it,kt) = (two_third * (-Q(3,i-1,j,k) + Q(3,i+1,j,k)) - one_twelfth * (-Q(3,i-2,j,k) + Q(3,i+2,j,k))) * dx(i)
      endif
      if (i <= nx .and. 1 <= j .and. j <= ny .and. 3 <= k .and. k <= nz-2) then
        vz(jj,it,kt) = (two_third * (-Q(3,i,j,k-1) + Q(3,i,j,k+1)) - one_twelfth * (-Q(3,i,j,k-2) + Q(3,i,j,k+2))) * dz(k)
        wz(jj,it,kt) = (two_third * (-Q(4,i,j,k-1) + Q(4,i,j,k+1)) - one_twelfth * (-Q(4,i,j,k-2) + Q(4,i,j,k+2))) * dz(k)
      endif
    enddo
    call syncthreads()
    j  = (blockIdx%y-1)*blockDim%y + jt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    if (3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-3 .and. 3 <= k .and. k <= nz-2) then
      block
        real(8), device :: mu3(3)
        mu3(:) = 0.0625d0 * (9.d0 * (mu(i,j-1:j+1,k) + mu(i,j:j+2,k)) - (mu(i,j-2:j,k) + mu(i,j+1:j+3,k)))
        block ! dQdy
          real(8), device :: kTy3(3)
          kTy3(:) = Cp_over_Pr * mu3(:) * &
                    (1.125d0 * (-T(i,j-1:j+1,k) + T(i,j:j+2,k)) - (-T(i,j-2:j,k) + T(i,j+1:j+3,k)) * one_24) * dy(j)
          kTy     = flux4(kTy3(:))
        end block
        call calc_tau_straight(mu3, v(jt-2:jt+3,it,kt), wz(jt-2:jt+3,it,kt), ux(jt-2:jt+3,it,kt), dy(j), tyy, vtyy)
        call calc_tau_cross(mu3, u(jt-2:jt+3,it,kt), vx(jt-2:jt+3,it,kt), dy(j), tyx, utyx)
        call calc_tau_cross(mu3, w(jt-2:jt+3,it,kt), vz(jt-2:jt+3,it,kt), dy(j), tyz, wtyz)
      end block
    else
      block
        real(8) my, muy, mvy, mwy, mvz, mwz, mux, mvx
        my  = 0.5d0 * (mu(i,j,k) + mu(i,j+1,k))
        kTy = Cp_over_Pr * my * (-T(i,j,k) + T(i,j+1,k)) * dy(j)
        block
          real(8), device :: mx(2)
          mx(1) = 0.25d0 * (mu(i-1,j,k) + mu(i,  j,k) + mu(i-1,j+1,k) + mu(i,  j+1,k))
          mx(2) = 0.25d0 * (mu(i,  j,k) + mu(i+1,j,k) + mu(i,  j+1,k) + mu(i+1,j+1,k))
          mux = 0.25d0 * (mx(1) * (-Q(2,i-1,j,k) + Q(2,i,j,k) - Q(2,i-1,j+1,k) + Q(2,i,j+1,k)) &
                        + mx(2) * (-Q(2,i,j,k) + Q(2,i+1,j,k) - Q(2,i,j+1,k) + Q(2,i+1,j+1,k))) * dx(i)
          mvx = 0.25d0 * (mx(1) * (-Q(3,i-1,j,k) + Q(3,i,j,k) - Q(3,i-1,j+1,k) + Q(3,i,j+1,k)) &
                        + mx(2) * (-Q(3,i,j,k) + Q(3,i+1,j,k) - Q(3,i,j+1,k) + Q(3,i+1,j+1,k))) * dx(i)
        end block
        my  = 0.5d0 * (mu(i,j,k) + mu(i,j+1,k))
        kTy = Cp_over_Pr * my * (-T(i,j,k) + T(i,j+1,k)) * dy(j)
        block
          real(8), device :: mz(2)
          mz(1) = 0.25d0 * (mu(i,j,k-1) + mu(i,j,k  ) + mu(i,j+1,k-1) + mu(i,j+1,k  ))
          mz(2) = 0.25d0 * (mu(i,j,k  ) + mu(i,j,k+1) + mu(i,j+1,k  ) + mu(i,j+1,k+1))
          mvz = 0.25d0 * (mz(1) * (-Q(3,i,j,k-1) + Q(3,i,j,k) - Q(3,i,j+1,k-1) + Q(3,i,j+1,k)) &
                        + mz(2) * (-Q(3,i,j,k) + Q(3,i,j,k+1) - Q(3,i,j+1,k) + Q(3,i,j+1,k+1))) * dz(k)
          mwz = 0.25d0 * (mz(1) * (-Q(4,i,j,k-1) + Q(4,i,j,k) - Q(4,i,j+1,k-1) + Q(4,i,j+1,k)) &
                        + mz(2) * (-Q(4,i,j,k) + Q(4,i,j,k+1) - Q(4,i,j+1,k) + Q(4,i,j+1,k+1))) * dz(k)
        end block
        muy  = my * (-u(jt,it,kt) + u(jt+1,it,kt)) * dy(j)
        mvy  = my * (-v(jt,it,kt) + v(jt+1,it,kt)) * dy(j)
        mwy  = my * (-w(jt,it,kt) + w(jt+1,it,kt)) * dy(j)
        tyx  = muy + mvx
        tyy  = two_third * (2.d0 * mvy - mwz - mux)
        tyz  = mvz + mwy
        utyx = 0.5d0 * (u(jt,it,kt) + u(jt+1,it,kt)) * tyx
        vtyy = 0.5d0 * (v(jt,it,kt) + v(jt+1,it,kt)) * tyy
        wtyz = 0.5d0 * (w(jt,it,kt) + w(jt+1,it,kt)) * tyz
      end block
    endif
    F(2,i-1,j,k-1) = F(2,i-1,j,k-1) - tyx
    F(3,i-1,j,k-1) = F(3,i-1,j,k-1) - tyy
    F(4,i-1,j,k-1) = F(4,i-1,j,k-1) - tyz
    F(5,i-1,j,k-1) = F(5,i-1,j,k-1) - (utyx + vtyy + wtyz + kTy)
  end subroutine calc_Fv4
 

  attributes(global) subroutine calc_Gv4(nx, ny, nz, dx, dy, dz, Q, T, mu, G)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz), T(nx,ny,nz), mu(nx,ny,nz)
    real(8), intent(inout), device :: G(5,nx-2,ny-2,nz-1)
    real(8), shared ::  u(-2:threadsGv%z+3,threadsGv%y,threadsGv%x)
    real(8), shared ::  v(-2:threadsGv%z+3,threadsGv%y,threadsGv%x)
    real(8), shared ::  w(-2:threadsGv%z+3,threadsGv%y,threadsGv%x)
    real(8), shared :: wx(-2:threadsGv%z+3,threadsGv%y,threadsGv%x)
    real(8), shared :: wy(-2:threadsGv%z+3,threadsGv%y,threadsGv%x)
    real(8), shared :: ux(-2:threadsGv%z+3,threadsGv%y,threadsGv%x)
    real(8), shared :: vy(-2:threadsGv%z+3,threadsGv%y,threadsGv%x)
    integer i, j, k, it, jt, kt, kk, k_base
    real(8) :: tzx, tzy, tzz, utzx, vtzy, wtzz, kTz
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k_base = (blockIdx%z-1)*blockDim%z
    do kk = kt-2, threadsGv%z+3, blockDim%z
      k = k_base + kk
      if (i <= nx .and. j <= ny .and. 1 <= k .and. k <= nz) then
        u(kk,jt,it) = Q(2,i,j,k)
        v(kk,jt,it) = Q(3,i,j,k)
        w(kk,jt,it) = Q(4,i,j,k)
      endif
      if (3 <= i .and. i <= nx-2 .and. j <= ny .and. 1 <= k .and. k <= nz) then
        ux(kk,jt,it) = (two_third * (-Q(2,i-1,j,k) + Q(2,i+1,j,k)) - one_twelfth * (-Q(2,i-2,j,k) + Q(2,i+2,j,k))) * dx(i)
        wx(kk,jt,it) = (two_third * (-Q(4,i-1,j,k) + Q(4,i+1,j,k)) - one_twelfth * (-Q(4,i-2,j,k) + Q(4,i+2,j,k))) * dx(i)
      endif
      if (i <= nx .and. 3 <= j .and. j <= ny-2 .and. 1 <= k .and. k <= nz) then
        vy(kk,jt,it) = (two_third * (-Q(3,i,j-1,k) + Q(3,i,j+1,k)) - one_twelfth * (-Q(3,i,j-2,k) + Q(3,i,j+2,k))) * dy(j)
        wy(kk,jt,it) = (two_third * (-Q(4,i,j-1,k) + Q(4,i,j+1,k)) - one_twelfth * (-Q(4,i,j-2,k) + Q(4,i,j+2,k))) * dy(j)
      endif
    enddo
    call syncthreads()
    k  = (blockIdx%z-1)*blockDim%z + kt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    if (3 <= i .and. i <= nx-2 .and. 3 <= j .and. j <= ny-2 .and. 3 <= k .and. k <= nz-3) then
      block
        real(8), device :: mu3(3)
        mu3(:) = 0.0625d0 * (9.d0 * (mu(i,j,k-1:k+1) + mu(i,j,k:k+2)) - (mu(i,j,k-2:k) + mu(i,j,k+1:k+3)))
        block ! dQdz
          real(8), device :: kTz3(3)
          kTz3(:) = Cp_over_Pr * mu3(:) * &
                    (1.125d0 * (-T(i,j,k-1:k+1) + T(i,j,k:k+2)) - (-T(i,j,k-2:k) + T(i,j,k+1:k+3)) * one_24) * dz(k)
          kTz     = flux4(kTz3(:))
        end block
        call calc_tau_straight(mu3, w(kt-2:kt+3,jt,it), ux(kt-2:kt+3,jt,it), vy(kt-2:kt+3,jt,it), dz(k), tzz, wtzz)
        call calc_tau_cross(mu3, u(kt-2:kt+3,jt,it), wx(kt-2:kt+3,jt,it), dz(k), tzx, utzx)
        call calc_tau_cross(mu3, v(kt-2:kt+3,jt,it), wy(kt-2:kt+3,jt,it), dz(k), tzy, vtzy)
      end block
    else
      block
        real(8) mz, muz, mvz, mwz, mwx, mux, mvy, mwy
        mz  = 0.5d0 * (mu(i,j,k) + mu(i,j,k+1))
        kTz = Cp_over_Pr * mz * (-T(i,j,k) + T(i,j,k+1)) * dz(k)
        block
          real(8), device :: mx(2)
          mx(1) = 0.25d0 * (mu(i-1,j,k) + mu(i,  j,k) + mu(i-1,j,k+1) + mu(i,  j,k+1))
          mx(2) = 0.25d0 * (mu(i,  j,k) + mu(i+1,j,k) + mu(i,  j,k+1) + mu(i+1,j,k+1))
          mux = 0.25d0 * (mx(1) * (-Q(2,i-1,j,k) + Q(2,i,j,k) - Q(2,i-1,j,k+1) + Q(2,i,j,k+1)) &
                        + mx(2) * (-Q(2,i,j,k) + Q(2,i+1,j,k) - Q(2,i,j,k+1) + Q(2,i+1,j,k+1))) * dx(i)
          mwx = 0.25d0 * (mx(1) * (-Q(4,i-1,j,k) + Q(4,i,j,k) - Q(4,i-1,j,k+1) + Q(4,i,j,k+1)) &
                        + mx(2) * (-Q(4,i,j,k) + Q(4,i+1,j,k) - Q(4,i,j,k+1) + Q(4,i+1,j,k+1))) * dx(i)
        end block
        block
          real(8), device :: my(2)
          my(1) = 0.25d0 * (mu(i,j-1,k) + mu(i,j,  k) + mu(i,j-1,k+1) + mu(i,j,  k+1))
          my(2) = 0.25d0 * (mu(i,j,  k) + mu(i,j+1,k) + mu(i,j,  k+1) + mu(i,j+1,k+1))
          mvy = 0.25d0 * (my(1) * (-Q(3,i,j-1,k) + Q(3,i,j,k) - Q(3,i,j-1,k+1) + Q(3,i,j,k+1)) &
                        + my(2) * (-Q(3,i,j,k) + Q(3,i,j+1,k) - Q(3,i,j,k+1) + Q(3,i,j+1,k+1))) * dy(j)
          mwy = 0.25d0 * (my(1) * (-Q(4,i,j-1,k) + Q(4,i,j,k) - Q(4,i,j-1,k+1) + Q(4,i,j,k+1)) &
                        + my(2) * (-Q(4,i,j,k) + Q(4,i,j+1,k) - Q(4,i,j,k+1) + Q(4,i,j+1,k+1))) * dy(j)
        end block
        mz   = 0.5d0 * (mu(i,j,k) + mu(i,j,k+1))
        kTz  = Cp_over_Pr * mz * (-T(i,j,k) + T(i,j,k+1)) * dz(k)
        muz  = mz * (-u(kt,jt,it) + u(kt+1,jt,it)) * dz(k)
        mvz  = mz * (-v(kt,jt,it) + v(kt+1,jt,it)) * dz(k)
        mwz  = mz * (-w(kt,jt,it) + w(kt+1,jt,it)) * dz(k)
        tzx  = mwx + muz
        tzy  = mvz + mwy
        tzz  = two_third * (2.d0 * mwz - mux - mvy)
        utzx = 0.5d0 * (u(kt,jt,it) + u(kt+1,jt,it)) * tzx
        vtzy = 0.5d0 * (v(kt,jt,it) + v(kt+1,jt,it)) * tzy
        wtzz = 0.5d0 * (w(kt,jt,it) + w(kt+1,jt,it)) * tzz
      end block
    endif
    G(2,i-1,j-1,k) = G(2,i-1,j-1,k) - tzx
    G(3,i-1,j-1,k) = G(3,i-1,j-1,k) - tzy
    G(4,i-1,j-1,k) = G(4,i-1,j-1,k) - tzz
    G(5,i-1,j-1,k) = G(5,i-1,j-1,k) - (utzx + vtzy + wtzz + kTz)
  end subroutine calc_Gv4
end module calc_visc4

