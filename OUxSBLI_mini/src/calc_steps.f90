module calc_steps
  use cudafor
  use cooperative_groups
  use mod_globals, only : dt
  use mod_constant, only : one_sixth
  implicit none
contains
  attributes(global) subroutine calc_R(nx, ny, nz, dx, dy, dz, E, F, G, R)
    integer, intent(in), value                               :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device             :: dx
    real(8), intent(in), dimension(ny-1), device             :: dy
    real(8), intent(in), dimension(nz-1), device             :: dz
    real(8), intent(in), dimension(5,nx-1,ny-2,nz-2), device :: E
    real(8), intent(in), dimension(5,nx-2,ny-1,nz-2), device :: F
    real(8), intent(in), dimension(5,nx-2,ny-2,nz-1), device :: G
    real(8), intent(out), dimension(5,nx-2,ny-2,nz-2),device :: R
    real(8) dtdydz, dtdzdx, dtdxdy, dx_next, E_curr, E_next
    integer i, j, k, l, lane
    integer(8) tmp_bits
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    lane     = iand(threadIdx%x - 1, 31)
    tmp_bits = __shfl_down_sync(z'ffffffff', transfer(dx(i), 0_8), 1)
    dx_next  = transfer(tmp_bits, 0.0_8)
    if (lane == 31 .or. i == nx-2) then
      dx_next = dx(i+1)
    endif
    dtdydz = dt * 0.25d0 * (dy(j) + dy(j+1)) * (dz(k) + dz(k+1))
    dtdzdx = dt * 0.25d0 * (dz(k) + dz(k+1)) * (dx(i) + dx_next)
    dtdxdy = dt * 0.25d0 * (dx(i) + dx_next) * (dy(j) + dy(j+1))
    do l = 1, 5
      E_curr   = E(l,i,j,k)
      tmp_bits = __shfl_down_sync(z'ffffffff', transfer(E_curr, 0_8), 1)
      E_next   = transfer(tmp_bits, 0.0_8)
      if (lane == 31 .or. i == nx-2) then
        E_next = E(l,i+1,j,k)
      endif
      R(l,i,j,k) = &
      &   dtdydz * (-E_curr     + E_next) &
      & + dtdzdx * (-F(l,i,j,k) + F(l,i,j+1,k)) &
      & + dtdxdy * (-G(l,i,j,k) + G(l,i,j,k+1))
    enddo
  end subroutine calc_R


  attributes(global) subroutine calc_step1(nx, ny, nz, coef, dx, dy, dz, E, F, G, Q, Q2)
    integer, intent(in), value                               :: nx, ny, nz
    real(8), intent(in), value                               :: coef
    real(8), intent(in), dimension(nx-1), device             :: dx
    real(8), intent(in), dimension(ny-1), device             :: dy
    real(8), intent(in), dimension(nz-1), device             :: dz
    real(8), intent(in), dimension(5,nx-1,ny-2,nz-2), device :: E
    real(8), intent(in), dimension(5,nx-2,ny-1,nz-2), device :: F
    real(8), intent(in), dimension(5,nx-2,ny-2,nz-1), device :: G
    real(8), intent(in), dimension(5,nx,ny,nz), device       :: Q
    real(8), intent(out), dimension(5,nx,ny,nz), device      :: Q2
    real(8) R, dtdydz, dtdzdx, dtdxdy, dx_next, E_curr, E_next
    integer i, j, k, l, lane
    integer(8) tmp_bits
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    lane     = iand(threadIdx%x - 1, 31)
    tmp_bits = __shfl_down_sync(z'ffffffff', transfer(dx(i), 0_8), 1)
    dx_next  = transfer(tmp_bits, 0.0_8)
    if (lane == 31 .or. i == nx-2) then
      dx_next = dx(i+1)
    endif
    dtdydz = dt * 0.25d0 * (dy(j) + dy(j+1)) * (dz(k) + dz(k+1))
    dtdzdx = dt * 0.25d0 * (dz(k) + dz(k+1)) * (dx(i) + dx_next)
    dtdxdy = dt * 0.25d0 * (dx(i) + dx_next) * (dy(j) + dy(j+1))
    do l = 1, 5
      E_curr   = E(l,i,j,k)
      tmp_bits = __shfl_down_sync(z'ffffffff', transfer(E_curr, 0_8), 1)
      E_next   = transfer(tmp_bits, 0.0_8)
      if (lane == 31 .or. i == nx-2) then
        E_next = E(l,i+1,j,k)
      endif
      R = dtdydz * (-E_curr     + E_next) &
      & + dtdzdx * (-F(l,i,j,k) + F(l,i,j+1,k)) &
      & + dtdxdy * (-G(l,i,j,k) + G(l,i,j,k+1))
      Q2(l,i+1,j+1,k+1) = Q(l,i+1,j+1,k+1) - coef * R
    enddo
  end subroutine calc_step1


  attributes(global) subroutine calc_step(nx, ny, nz, coef1, coef2, dx, dy, dz, E, F, G, Q, Q2, Rs)
    integer, intent(in), value                                  :: nx, ny, nz
    real(8), intent(in), value                                  :: coef1, coef2
    real(8), intent(in), dimension(nx-1), device                :: dx
    real(8), intent(in), dimension(ny-1), device                :: dy
    real(8), intent(in), dimension(nz-1), device                :: dz
    real(8), intent(in), dimension(5,nx-1,ny-2,nz-2), device    :: E
    real(8), intent(in), dimension(5,nx-2,ny-1,nz-2), device    :: F
    real(8), intent(in), dimension(5,nx-2,ny-2,nz-1), device    :: G
    real(8), intent(in), dimension(5,nx,ny,nz), device          :: Q
    real(8), intent(out), dimension(5,nx,ny,nz), device         :: Q2
    real(8), intent(inout), dimension(5,nx-2,ny-2,nz-2), device :: Rs
    real(8) R, dtdydz, dtdzdx, dtdxdy, dx_next, E_curr, E_next
    integer i, j, k, l, lane
    integer(8) tmp_bits
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    lane     = iand(threadIdx%x - 1, 31)
    tmp_bits = __shfl_down_sync(z'ffffffff', transfer(dx(i), 0_8), 1)
    dx_next  = transfer(tmp_bits, 0.0_8)
    if (lane == 31 .or. i == nx-2) then
      dx_next = dx(i+1)
    endif
    dtdydz = dt * 0.25d0 * (dy(j) + dy(j+1)) * (dz(k) + dz(k+1))
    dtdzdx = dt * 0.25d0 * (dz(k) + dz(k+1)) * (dx(i) + dx_next)
    dtdxdy = dt * 0.25d0 * (dx(i) + dx_next) * (dy(j) + dy(j+1))
    do l = 1, 5
      E_curr   = E(l,i,j,k)
      tmp_bits = __shfl_down_sync(z'ffffffff', transfer(E_curr, 0_8), 1)
      E_next   = transfer(tmp_bits, 0.0_8)
      if (lane == 31 .or. i == nx-2) then
        E_next = E(l,i+1,j,k)
      endif
      R = dtdydz * (-E_curr     + E_next) &
      & + dtdzdx * (-F(l,i,j,k) + F(l,i,j+1,k)) &
      & + dtdxdy * (-G(l,i,j,k) + G(l,i,j,k+1))
      Q2(l,i+1,j+1,k+1) = Q(l,i+1,j+1,k+1) - coef1 * R
      Rs(l,i,j,k) = Rs(l,i,j,k) + coef2 * R
    enddo
  end subroutine calc_step
 

  attributes(global) subroutine calc_step2_3(nx, ny, nz, coef1, coef2, coef3, coef4, dx, dy, dz, E, F, G, Qin, Qout)
    integer, intent(in), value                               :: nx, ny, nz
    real(8), intent(in), value                               :: coef1, coef2, coef3, coef4
    real(8), intent(in), dimension(nx-1), device             :: dx
    real(8), intent(in), dimension(ny-1), device             :: dy
    real(8), intent(in), dimension(nz-1), device             :: dz
    real(8), intent(in), dimension(5,nx-1,ny-2,nz-2), device :: E
    real(8), intent(in), dimension(5,nx-2,ny-1,nz-2), device :: F
    real(8), intent(in), dimension(5,nx-2,ny-2,nz-1), device :: G
    real(8), intent(in), dimension(5,nx,ny,nz), device       :: Qin
    real(8), intent(inout), dimension(5,nx,ny,nz), device    :: Qout
    real(8) R, dtdydz, dtdzdx, dtdxdy, dx_next, E_curr, E_next
    integer i, j, k, l, lane
    integer(8) tmp_bits
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    lane     = iand(threadIdx%x - 1, 31)
    tmp_bits = __shfl_down_sync(z'ffffffff', transfer(dx(i), 0_8), 1)
    dx_next  = transfer(tmp_bits, 0.0_8)
    if (lane == 31 .or. i == nx-2) then
      dx_next = dx(i+1)
    endif
    dtdydz = dt * 0.25d0 * (dy(j) + dy(j+1)) * (dz(k) + dz(k+1))
    dtdzdx = dt * 0.25d0 * (dz(k) + dz(k+1)) * (dx(i) + dx_next)
    dtdxdy = dt * 0.25d0 * (dx(i) + dx_next) * (dy(j) + dy(j+1))
    do l = 1, 5
      E_curr   = E(l,i,j,k)
      tmp_bits = __shfl_down_sync(z'ffffffff', transfer(E_curr, 0_8), 1)
      E_next   = transfer(tmp_bits, 0.0_8)
      if (lane == 31 .or. i == nx-2) then
        E_next = E(l,i+1,j,k)
      endif
      R = dtdydz * (-E_curr     + E_next) &
      & + dtdzdx * (-F(l,i,j,k) + F(l,i,j+1,k)) &
      & + dtdxdy * (-G(l,i,j,k) + G(l,i,j,k+1))
      Qout(l,i+1,j+1,k+1) = (coef1 * Qin(l,i+1,j+1,k+1) + coef2 * Qout(l,i+1,j+1,k+1) - coef3 * R) / coef4
    enddo
  end subroutine calc_step2_3
  
 
  attributes(global) subroutine calc_step4(nx, ny, nz, dx, dy, dz, E, F, G, Rs, Q)
    integer, intent(in), value                                  :: nx, ny, nz
    real(8), intent(in), dimension(nx-1), device                :: dx
    real(8), intent(in), dimension(ny-1), device                :: dy
    real(8), intent(in), dimension(nz-1), device                :: dz
    real(8), intent(in), dimension(5,nx-1,ny-2,nz-2), device    :: E
    real(8), intent(in), dimension(5,nx-2,ny-1,nz-2), device    :: F
    real(8), intent(in), dimension(5,nx-2,ny-2,nz-1), device    :: G
    real(8), intent(inout), dimension(5,nx-2,ny-2,nz-2), device :: Rs
    real(8), intent(inout), dimension(5,nx,ny,nz), device       :: Q
    real(8) R, dtdydz, dtdzdx, dtdxdy, dx_next, E_curr, E_next
    integer i, j, k, l, lane
    integer(8) tmp_bits
    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
    if (nx-2 < i .or. ny-2 < j .or. nz-2 < k) return
    lane     = iand(threadIdx%x - 1, 31)
    tmp_bits = __shfl_down_sync(z'ffffffff', transfer(dx(i), 0_8), 1)
    dx_next  = transfer(tmp_bits, 0.0_8)
    if (lane == 31 .or. i == nx-2) then
      dx_next = dx(i+1)
    endif
    dtdydz = dt * 0.25d0 * (dy(j) + dy(j+1)) * (dz(k) + dz(k+1))
    dtdzdx = dt * 0.25d0 * (dz(k) + dz(k+1)) * (dx(i) + dx_next)
    dtdxdy = dt * 0.25d0 * (dx(i) + dx_next) * (dy(j) + dy(j+1))
    do l = 1, 5
      E_curr   = E(l,i,j,k)
      tmp_bits = __shfl_down_sync(z'ffffffff', transfer(E_curr, 0_8), 1)
      E_next   = transfer(tmp_bits, 0.0_8)
      if (lane == 31 .or. i == nx-2) then
        E_next = E(l,i+1,j,k)
      endif
      R = dtdydz * (-E_curr     + E_next) &
      & + dtdzdx * (-F(l,i,j,k) + F(l,i,j+1,k)) &
      & + dtdxdy * (-G(l,i,j,k) + G(l,i,j,k+1))
      Rs(l,i,j,k) = Rs(l,i,j,k) + R
      Q(l,i+1,j+1,k+1) = Q(l,i+1,j+1,k+1) - Rs(l,i,j,k) * one_sixth
      Rs(l,i,j,k) = 0.d0
    enddo
  end subroutine calc_step4
end module calc_steps

