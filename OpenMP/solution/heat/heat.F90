program heat

  implicit none 

  ! precision parameter for 64bit numbers
  INTEGER, PARAMETER :: DP = 8

  integer :: i,k,it 

  !     integer imax,kmax,itmax 
  !     parameter (imax=20,kmax=11)
  !     parameter (itmax=20000)

#ifndef imax
#  define imax  20
#endif 
#ifndef kmax
#  define kmax  20
#endif 
#ifndef itmax 
#  define itmax  15000
#endif 

#ifndef prtrows
#  define prtrows 10
#endif 
#ifndef prtcols
#  define prtcols 10
#endif 

  REAL (KIND=DP), parameter :: eps = 10.0**-08
  REAL (KIND=DP), dimension (0:imax,0:kmax) :: phi
  REAL (KIND=DP), dimension (1:imax-1,1:kmax-1) :: phin

  REAL (KIND=DP) :: dx,dy,dx2,dy2,dx2i,dy2i,dt,dphi,dphimax

  ! times using cpu_time
  real :: t0
  real :: t1

  !--unused-- include 'omp_lib.h' 
  !$    integer omp_get_num_threads 
  !$    REAL (KIND=DP) omp_get_wtime
  !$    REAL (KIND=DP) wt0,wt1   

  ! END DECLARATION


  !$omp parallel
  !$omp single
  !$    write(*,*)'OpenMP-parallel with',omp_get_num_threads(),'threads'
  !$omp end single
  !$omp end parallel

  dx = 1.0_DP / kmax
  dy = 1.0_DP / imax
  dx2 = dx*dx
  dy2 = dy*dy
  dx2i = 1.0_DP / dx2
  dy2i = 1.0_DP / dy2
  dt = min(dx2,dy2) / 4.0_DP

  ! start values 0.d0
  do i = 1, imax-1
     do k = 0, kmax-1
        phi(i,k) = 0.0_DP
     enddo
  enddo

  ! start values 1.d0
  do i=0, imax
     phi(i,kmax) = 1.0_DP
  enddo
!!$  phi(:,kmax) = 1.0_DP

  ! start values dx
  phi(0,0) = 0.0_DP
  phi(imax,0) = 0.0_DP

  do k = 1,kmax-1
     phi(0,k) = phi(0,k-1) + dx
     phi(imax,k) = phi(imax,k-1) + dx
  enddo

  ! print array
  write (*,'(/,a)') 'Heat Conduction 2d'
  write (*,'(/,4(a,1pg12.4))') 'dx =',dx,', dy =',dy,', dt =',dt,', eps =',eps
  write (*,'(/,a)') 'start values'

  call heatpr(phi,imax,kmax)

  !
  !$    wt0=omp_get_wtime() 

  call cpu_time(t0)

  ! iteration
  do it = 1,itmax
     dphimax = 0.0

     do i = 1,imax-1
        do k = 1,kmax-1
           dphi = (phi(i+1,k)+phi(i-1,k)-2.*phi(i,k))*dy2i &
                +(phi(i,k+1)+phi(i,k-1)-2.*phi(i,k))*dx2i
           dphi = dphi*dt
           dphimax = max(dphimax,dphi)
           phin(i,k) = phi(i,k)+dphi
        enddo
     enddo

     ! save values
     do i =1,imax-1
        do k = 1,kmax-1
           phi(i,k) =phin(i,k)
        enddo
     enddo

     if(dphimax < eps) EXIT

  enddo

  !
  call cpu_time(t1)

  !$    wt1 = omp_get_wtime() 

  ! print array
  write (*,'(/,a,i6,a)') 'phi after',it,' iterations'

  call heatpr(phi,imax,kmax)

  write (*,'(/,a,1pg12.4)')  'cpu_time  :   ', t1-t0
  !$    write (*,'(/,a,1pg12.4)')  'omp_get_wtime:', wt1-wt0


  stop

end program heat


!
!
!
!
subroutine heatpr(phi,imx,kmx)

  implicit none

  ! precision parameter for 64bit numbers
  INTEGER, PARAMETER :: DP = 8

  REAL (KIND=DP), dimension (0:imx,0:kmx) :: phi

  integer :: imx, kmx 
  integer :: i, prt_i,prt_k, prt_i_max, prt_k_max
  integer :: k(20) 


  prt_i_max = min(prtrows, imx+1) 
  prt_k_max = min(prtcols, kmx+1, 20) 

  do prt_k = 0, prt_k_max-1
     k(prt_k) = 1.0*prt_k/(prt_k_max-1)*kmx + 0.5/*rounding*/
  end do

  write (*,'(''i  |  k='',20(''  '',I6))') (k(prt_k), prt_k = 0,prt_k_max-1)
  write (*,'(''-------+'',20A)') ('--------', prt_k = 0,prt_k_max-1)

  do prt_i = 0, prt_i_max-1
     i = 1.0*prt_i/(prt_i_max-1)*imx + 0.5/*rounding*/
     write (*, '(I6,'' |'',20(2X,F6.4))') &
          i, (phi(i,k(prt_k)), prt_k = 0,prt_k_max-1) 
  end do


  return

end subroutine heatpr
