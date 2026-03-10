module calc_visc2
  use mod_globals, only : id_LL, gamma, R, Pr, Prt, dt, threadsEv, threadsFv, threadsGv
  use mod_constant, only : Cp, gamma_1, Cp_over_Pr, one_third, two_third
  implicit none
contains
  attributes(global) subroutine calc_Ev2(nx, ny, nz, dx, dy, dz, Q, T, mu, E)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz), T(nx,ny,nz), mu(nx,ny,nz)
    real(8), intent(inout), device :: E(5,nx-1,ny-2,nz-2)
    real(8), shared :: u(threadsEv%x+1,0:threadsEv%y+1,0:threadsEv%z+1)
    real(8), shared :: v(threadsEv%x+1,0:threadsEv%y+1,threadsEv%z)
    real(8), shared :: w(threadsEv%x+1,0:threadsEv%z+1,threadsEv%y)
    integer i, j, k, it, jt, kt
    real(8) :: txx, txy, txz, utxx, vtxy, wtxz, kTx
    real(8) m1, m2, mx, mux, mvx, mwx, muy, mvy, muz, mwz
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    u(it,jt-1:jt+1,kt-1:kt+1) = Q(2,i,j-1:j+1,k-1:k+1)
    v(it,jt-1:jt+1,kt)        = Q(3,i,j-1:j+1,k)
    w(it,kt-1:kt+1,jt)        = Q(4,i,j,k-1:k+1)
    if (it == blockDim%x) then
      u(it+1,jt-1:jt+1,kt-1:kt+1) = Q(2,i+1,j-1:j+1,k-1:k+1)
      v(it+1,jt-1:jt+1,kt)        = Q(3,i+1,j-1:j+1,k)
      w(it+1,kt-1:kt+1,jt)        = Q(4,i+1,j,k-1:k+1)
    endif 
    call syncthreads()

    mx  = 0.5d0 * (mu(i,j,k) + mu(i+1,j,k))
    kTx = Cp_over_Pr * mx * (-T(i,j,k) + T(i+1,j,k)) * dx(i)
    block
      real(8) my1, my2
      my1 = 0.25d0 * (mu(i,j-1,k) + mu(i,j,  k) + mu(i+1,j-1,k) + mu(i+1,j,  k))
      my2 = 0.25d0 * (mu(i,j,  k) + mu(i,j+1,k) + mu(i+1,j,  k) + mu(i+1,j+1,k))
      muy = 0.25d0 * (my1 * (-u(it,jt-1,kt) - u(it+1,jt-1,kt)) + (my1 - my2) * (u(it,jt,kt) + u(it+1,jt,kt)) &
                    + my2 * ( u(it,jt+1,kt) + u(it+1,jt+1,kt))) * dy(j)
      mvy = 0.25d0 * (my1 * (-v(it,jt-1,kt) - v(it+1,jt-1,kt)) + (my1 - my2) * (v(it,jt,kt) + v(it+1,jt,kt)) &
                    + my2 * ( v(it,jt+1,kt) + v(it+1,jt+1,kt))) * dy(j)
    end block
    block
      real(8) mz1, mz2
      mz1 = 0.25d0 * (mu(i,j,k-1) + mu(i,j,k  ) + mu(i+1,j,k-1) + mu(i+1,j,k  ))
      mz2 = 0.25d0 * (mu(i,j,k  ) + mu(i,j,k+1) + mu(i+1,j,k  ) + mu(i+1,j,k+1))
      muz = 0.25d0 * (mz1 * (-u(it,jt,kt-1) - u(it+1,jt,kt-1)) + (mz1 - mz2) * (u(it,jt,kt) + u(it+1,jt,kt)) &
                    + mz2 * ( u(it,jt,kt+1) + u(it+1,jt,kt+1))) * dz(k)
      mwz = 0.25d0 * (mz1 * (-w(it,kt-1,jt) - w(it+1,kt-1,jt)) + (mz1 - mz2) * (w(it,kt,jt) + w(it+1,kt,jt)) &
                    + mz2 * ( w(it,kt+1,jt) + w(it+1,kt+1,jt))) * dz(k)
    end block
    mux = mx * (-u(it,jt,kt) + u(it+1,jt,kt)) * dx(i)
    mvx = mx * (-v(it,jt,kt) + v(it+1,jt,kt)) * dx(i)
    mwx = mx * (-w(it,kt,jt) + w(it+1,kt,jt)) * dx(i)
    txx = two_third * (2.d0 * mux - mvy - mwz)
    txy = muy + mvx
    txz = mwx + muz
    utxx = 0.5d0 * (u(it,jt,kt) + u(it+1,jt,kt)) * txx
    vtxy = 0.5d0 * (v(it,jt,kt) + v(it+1,jt,kt)) * txy
    wtxz = 0.5d0 * (w(it,kt,jt) + w(it+1,kt,jt)) * txz
    E(2,i,j-1,k-1) = E(2,i,j-1,k-1) - txx
    E(3,i,j-1,k-1) = E(3,i,j-1,k-1) - txy
    E(4,i,j-1,k-1) = E(4,i,j-1,k-1) - txz
    E(5,i,j-1,k-1) = E(5,i,j-1,k-1) - (utxx + vtxy + wtxz + kTx)
  end subroutine calc_Ev2
 

  attributes(global) subroutine calc_Fv2(nx, ny, nz, dy, dx, dz, Q, T, mu, F)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz), T(nx,ny,nz), mu(nx,ny,nz)
    real(8), intent(inout), device :: F(5,nx-2,ny-1,nz-2)
    real(8), shared :: u(threadsFv%y+1,0:threadsFv%x+1,threadsFv%z)
    real(8), shared :: v(threadsFv%y+1,0:threadsFv%x+1,0:threadsFv%z+1)
    real(8), shared :: w(threadsFv%y+1,0:threadsFv%z+1,threadsFv%x)
    integer i, j, k, it, jt, kt
    real(8) :: tyx, tyy, tyz, utyx, vtyy, wtyz, kTy
    real(8) m1, m2, my, muy, mvy, mwy, mvz, mwz, mux, mvx
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt
    k  = (blockIdx%z-1)*blockDim%z + kt + 1
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    u(jt,it-1:it+1,kt)        = Q(2,i-1:i+1,j,k)
    v(jt,it-1:it+1,kt-1:kt+1) = Q(3,i-1:i+1,j,k-1:k+1)
    w(jt,kt-1:kt+1,it)        = Q(4,i,j,k-1:k+1)
    if (jt == blockDim%y) then
      u(jt+1,it-1:it+1,kt)        = Q(2,i-1:i+1,j+1,k)
      v(jt+1,it-1:it+1,kt-1:kt+1) = Q(3,i-1:i+1,j+1,k-1:k+1)
      w(jt+1,kt-1:kt+1,it)        = Q(4,i,j+1,k-1:k+1)
    endif 
    call syncthreads()

    block
      real(8) mx1, mx2
      mx1 = 0.25d0 * (mu(i-1,j,k) + mu(i,  j,k) + mu(i-1,j+1,k) + mu(i,  j+1,k))
      mx2 = 0.25d0 * (mu(i,  j,k) + mu(i+1,j,k) + mu(i,  j+1,k) + mu(i+1,j+1,k))
      mux = 0.25d0 * (mx1 * (-u(jt,it-1,kt) - u(jt+1,it-1,kt)) + (mx1 - mx2) * (u(jt,it,kt) + u(jt+1,it,kt)) &
                    + mx2 * ( u(jt,it+1,kt) + u(jt+1,it+1,kt))) * dx(i)
      mvx = 0.25d0 * (mx1 * (-v(jt,it-1,kt) - v(jt+1,it-1,kt)) + (mx1 - mx2) * (v(jt,it,kt) + v(jt+1,it,kt)) &
                    + mx2 * ( v(jt,it+1,kt) + v(jt+1,it+1,kt))) * dx(i)
    end block
    my  = 0.5d0 * (mu(i,j,k) + mu(i,j+1,k))
    kTy = Cp_over_Pr * my * (-T(i,j,k) + T(i,j+1,k)) * dy(j)
    block
      real(8) mz1, mz2
      mz1 = 0.25d0 * (mu(i,j,k-1) + mu(i,j,k  ) + mu(i,j+1,k-1) + mu(i,j+1,k  ))
      mz2 = 0.25d0 * (mu(i,j,k  ) + mu(i,j,k+1) + mu(i,j+1,k  ) + mu(i,j+1,k+1))
      mvz = 0.25d0 * (mz1 * (-v(jt,it,kt-1) - v(jt+1,it,kt-1)) + (mz1 - mz2) * (v(jt,it,kt) + v(jt+1,it,kt)) &
                    + mz2 * ( v(jt,it,kt+1) + v(jt+1,it,kt+1))) * dz(k)
      mwz = 0.25d0 * (mz1 * (-w(jt,kt-1,it) - w(jt+1,kt-1,it)) + (mz1 - mz2) * (w(jt,kt,it) + w(jt+1,kt,it)) &
                    + mz2 * ( w(jt,kt+1,it) + w(jt+1,kt+1,it))) * dz(k)
    end block
    muy = my * (-u(jt,it,kt) + u(jt+1,it,kt)) * dy(j)
    mvy = my * (-v(jt,it,kt) + v(jt+1,it,kt)) * dy(j)
    mwy = my * (-w(jt,kt,it) + w(jt+1,kt,it)) * dy(j)
    tyx = muy + mvx
    tyy = two_third * (2.d0 * mvy - mwz - mux)
    tyz = mvz + mwy
    utyx = 0.5d0 * (u(jt,it,kt) + u(jt+1,it,kt)) * tyx
    vtyy = 0.5d0 * (v(jt,it,kt) + v(jt+1,it,kt)) * tyy
    wtyz = 0.5d0 * (w(jt,kt,it) + w(jt+1,kt,it)) * tyz
    F(2,i-1,j,k-1) = F(2,i-1,j,k-1) - tyx
    F(3,i-1,j,k-1) = F(3,i-1,j,k-1) - tyy
    F(4,i-1,j,k-1) = F(4,i-1,j,k-1) - tyz
    F(5,i-1,j,k-1) = F(5,i-1,j,k-1) - (utyx + vtyy + wtyz + kTy)
  end subroutine calc_Fv2


  attributes(global) subroutine calc_Gv2(nx, ny, nz, dx, dy, dz, Q, T, mu, G)
    integer, intent(in), value     :: nx, ny, nz
    real(8), intent(in), device    :: dx(nx-1) ! 1 / dx
    real(8), intent(in), device    :: dy(ny-1) ! 1 / dy
    real(8), intent(in), device    :: dz(nz-1) ! 1 / dz
    real(8), intent(in), device    :: Q(5,nx,ny,nz), T(nx,ny,nz), mu(nx,ny,nz)
    real(8), intent(inout), device :: G(5,nx-2,ny-2,nz-1)
    real(8), shared :: u(threadsGv%z+1,0:threadsGv%x+1,threadsGv%y)
    real(8), shared :: v(threadsGv%z+1,0:threadsGv%y+1,threadsGv%x)
    real(8), shared :: w(threadsGv%z+1,0:threadsGv%x+1,0:threadsGv%y+1)
    integer i, j, k, it, jt, kt
    real(8) :: tzx, tzy, tzz, utzx, vtzy, wtzz, kTz
    real(8) m1, m2, mz, muz, mvz, mwz, mwx, mux, mvy, mwy
    it = threadIdx%x
    jt = threadIdx%y
    kt = threadIdx%z
    i  = (blockIdx%x-1)*blockDim%x + it + 1
    j  = (blockIdx%y-1)*blockDim%y + jt + 1
    k  = (blockIdx%z-1)*blockDim%z + kt
    if (nx-1 < i .or. ny-1 < j .or. nz-1 < k) return
    u(kt,it-1:it+1,jt)        = Q(2,i-1:i+1,j,k)
    v(kt,jt-1:jt+1,it)        = Q(3,i,j-1:j+1,k)
    w(kt,it-1:it+1,jt-1:jt+1) = Q(4,i-1:i+1,j-1:j+1,k)
    if (kt == blockDim%z) then
      u(kt+1,it-1:it+1,jt)        = Q(2,i-1:i+1,j,k+1)
      v(kt+1,jt-1:jt+1,it)        = Q(3,i,j-1:j+1,k+1)
      w(kt+1,it-1:it+1,jt-1:jt+1) = Q(4,i-1:i+1,j-1:j+1,k+1)
    endif 
    call syncthreads()

    block
      real(8) mx1, mx2
      mx1 = 0.25d0 * (mu(i-1,j,k) + mu(i,  j,k) + mu(i-1,j,k+1) + mu(i,  j,k+1))
      mx2 = 0.25d0 * (mu(i,  j,k) + mu(i+1,j,k) + mu(i,  j,k+1) + mu(i+1,j,k+1))
      mux = 0.25d0 * (mx1 * (-u(kt,it-1,jt) - u(kt+1,it-1,jt)) + (mx1 - mx2) * (u(kt,it,jt) + u(kt+1,it,jt)) &
                    + mx2 * ( u(kt,it+1,jt) + u(kt+1,it+1,jt))) * dx(i)
      mwx = 0.25d0 * (mx1 * (-w(kt,it-1,jt) - w(kt+1,it-1,jt)) + (mx1 - mx2) * (w(kt,it,jt) + w(kt+1,it,jt)) &
                    + mx2 * ( w(kt,it+1,jt) + w(kt+1,it+1,jt))) * dx(i)
    end block
    block
      real(8) my1, my2
      my1 = 0.25d0 * (mu(i,j-1,k) + mu(i,j,  k) + mu(i,j-1,k+1) + mu(i,j,  k+1))
      my2 = 0.25d0 * (mu(i,j,  k) + mu(i,j+1,k) + mu(i,j,  k+1) + mu(i,j+1,k+1))
      mvy = 0.25d0 * (my1 * (-v(kt,jt-1,it) - v(kt+1,jt-1,it)) + (my1 - my2) * (v(kt,jt,it) + v(kt+1,jt,it)) &
                    + my2 * ( v(kt,jt+1,it) + v(kt+1,jt+1,it))) * dy(j)
      mwy = 0.25d0 * (my1 * (-w(kt,it,jt-1) - w(kt+1,it,jt-1)) + (my1 - my2) * (w(kt,it,jt) + w(kt+1,it,jt)) &
                    + my2 * ( w(kt,it,jt+1) + w(kt+1,it,jt+1))) * dy(j)
    end block
    mz  = 0.5d0 * (mu(i,j,k) + mu(i,j,k+1))
    kTz = Cp_over_Pr * mz * (-T(i,j,k) + T(i,j,k+1)) * dz(k)
    muz = mz * (-u(kt,it,jt) + u(kt+1,it,jt)) * dz(k)
    mvz = mz * (-v(kt,jt,it) + v(kt+1,jt,it)) * dz(k)
    mwz = mz * (-w(kt,it,jt) + w(kt+1,it,jt)) * dz(k)
    tzx = mwx + muz
    tzy = mvz + mwy
    tzz = two_third * (2.d0 * mwz - mux - mvy)
    utzx = 0.5d0 * (u(kt,it,jt) + u(kt+1,it,jt)) * tzx
    vtzy = 0.5d0 * (v(kt,jt,it) + v(kt+1,jt,it)) * tzy
    wtzz = 0.5d0 * (w(kt,it,jt) + w(kt+1,it,jt)) * tzz
    G(2,i-1,j-1,k) = G(2,i-1,j-1,k) - tzx
    G(3,i-1,j-1,k) = G(3,i-1,j-1,k) - tzy
    G(4,i-1,j-1,k) = G(4,i-1,j-1,k) - tzz
    G(5,i-1,j-1,k) = G(5,i-1,j-1,k) - (utzx + vtzy + wtzz + kTz)
  end subroutine calc_Gv2
end module calc_visc2

