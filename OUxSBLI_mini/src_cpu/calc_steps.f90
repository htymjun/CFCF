module calc_steps
  use mod_globals, only : dt
  use mod_constant, only : one_sixth
  implicit none
contains
  subroutine calc_step(nx, ny, nz, coef1, coef2, dx, dy, dz, E, F, G, Q, Q2, Rs)
    integer, intent(in)                                  :: nx, ny, nz
    real(8), intent(in)                                  :: coef1, coef2
    real(8), intent(in), dimension(nx-1)                 :: dx
    real(8), intent(in), dimension(ny-1)                 :: dy
    real(8), intent(in), dimension(nz-1)                 :: dz
    real(8), intent(in), dimension(5,nx-1,ny-2,nz-2)     :: E
    real(8), intent(in), dimension(5,nx-2,ny-1,nz-2)     :: F
    real(8), intent(in), dimension(5,nx-2,ny-2,nz-1)     :: G
    real(8), intent(in), dimension(5,nx,ny,nz)            :: Q
    real(8), intent(out), dimension(5,nx,ny,nz)           :: Q2
    real(8), intent(inout), dimension(5,nx-2,ny-2,nz-2)  :: Rs
    real(8) Rv, dtdydz, dtdzdx, dtdxdy
    integer i, j, k, l
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j,k,l,dtdydz,dtdzdx,dtdxdy,Rv)
    do k = 1, nz-2
      do j = 1, ny-2
        do i = 1, nx-2
          dtdydz = dt * 0.25d0 * (dy(j) + dy(j+1)) * (dz(k) + dz(k+1))
          dtdzdx = dt * 0.25d0 * (dz(k) + dz(k+1)) * (dx(i) + dx(i+1))
          dtdxdy = dt * 0.25d0 * (dx(i) + dx(i+1)) * (dy(j) + dy(j+1))
          do l = 1, 5
            Rv = dtdydz * (-E(l,i,j,k) + E(l,i+1,j,k)) &
            &  + dtdzdx * (-F(l,i,j,k) + F(l,i,j+1,k)) &
            &  + dtdxdy * (-G(l,i,j,k) + G(l,i,j,k+1))
            Q2(l,i+1,j+1,k+1) = Q(l,i+1,j+1,k+1) - coef1 * Rv
            Rs(l,i,j,k) = Rs(l,i,j,k) + coef2 * Rv
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine calc_step


  subroutine calc_step4(nx, ny, nz, dx, dy, dz, E, F, G, Rs, Q)
    integer, intent(in)                                  :: nx, ny, nz
    real(8), intent(in), dimension(nx-1)                 :: dx
    real(8), intent(in), dimension(ny-1)                 :: dy
    real(8), intent(in), dimension(nz-1)                 :: dz
    real(8), intent(in), dimension(5,nx-1,ny-2,nz-2)     :: E
    real(8), intent(in), dimension(5,nx-2,ny-1,nz-2)     :: F
    real(8), intent(in), dimension(5,nx-2,ny-2,nz-1)     :: G
    real(8), intent(inout), dimension(5,nx-2,ny-2,nz-2)  :: Rs
    real(8), intent(inout), dimension(5,nx,ny,nz)         :: Q
    real(8) Rv, dtdydz, dtdzdx, dtdxdy
    integer i, j, k, l
    !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j,k,l,dtdydz,dtdzdx,dtdxdy,Rv)
    do k = 1, nz-2
      do j = 1, ny-2
        do i = 1, nx-2
          dtdydz = dt * 0.25d0 * (dy(j) + dy(j+1)) * (dz(k) + dz(k+1))
          dtdzdx = dt * 0.25d0 * (dz(k) + dz(k+1)) * (dx(i) + dx(i+1))
          dtdxdy = dt * 0.25d0 * (dx(i) + dx(i+1)) * (dy(j) + dy(j+1))
          do l = 1, 5
            Rv = dtdydz * (-E(l,i,j,k) + E(l,i+1,j,k)) &
            &  + dtdzdx * (-F(l,i,j,k) + F(l,i,j+1,k)) &
            &  + dtdxdy * (-G(l,i,j,k) + G(l,i,j,k+1))
            Rs(l,i,j,k) = Rs(l,i,j,k) + Rv
            Q(l,i+1,j+1,k+1) = Q(l,i+1,j+1,k+1) - Rs(l,i,j,k) * one_sixth
            Rs(l,i,j,k) = 0.d0
          enddo
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine calc_step4
end module calc_steps
