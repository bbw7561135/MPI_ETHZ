! using critical section inside of main loop
! 
      program heat
      implicit none 
      integer i,k,it 
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
      double precision eps
      parameter (eps=1.d-08)
      double precision phi(0:imax,0:kmax), phin(1:imax-1,1:kmax-1)
      double precision dx,dy,dx2,dy2,dx2i,dy2i,dt,dphi,dphimax
! times using cpu_time
      real t0
      real t1
!--unused-- include 'omp_lib.h' 
!$    integer omp_get_num_threads 
!$    double precision omp_get_wtime
!$    double precision wt0,wt1   
!
!$omp parallel
!$omp single
!$    write(*,*)'OpenMP-parallel with',omp_get_num_threads(),'threads'
!$omp end single
!$omp end parallel
!
      dx=1.d0/kmax
      dy=1.d0/imax
      dx2=dx*dx
      dy2=dy*dy
      dx2i=1.d0/dx2
      dy2i=1.d0/dy2
      dt=min(dx2,dy2)/4.d0
! start values 0.d0
!$OMP PARALLEL PRIVATE(i,k), SHARED(phi)
!     'i'-loop should be inner loop because phi(i,k) is contiguous with index i
!$OMP DO
      do k=0,kmax-1
        do i=1,imax-1
          phi(i,k)=0.d0
        enddo
      enddo
!$OMP END DO
! start values 1.d0
!$OMP DO
      do i=0,imax
        phi(i,kmax)=1.d0
      enddo
!$OMP END DO
!$OMP END PARALLEL
! start values dx
      phi(0,0)=0.d0
      phi(imax,0)=0.d0
      do k=1,kmax-1
        phi(0,k)=phi(0,k-1)+dx
        phi(imax,k)=phi(imax,k-1)+dx
      enddo
! print array
      write (*,'(/,a)')
     f      'Heat Conduction 2d'
      write (*,'(/,4(a,1pg12.4))')
     f      'dx =',dx,', dy =',dy,', dt =',dt,', eps =',eps
      write (*,'(/,a)')
     f      'start values'
      call heatpr(phi,imax,kmax)
!
!$    wt0=omp_get_wtime() 
      call cpu_time(t0)
! iteration
      do it=1,itmax
       dphimax=0.
!$OMP  PARALLEL PRIVATE(i,k,dphi)
!$OMP&          SHARED(phi,phin,dx2i,dy2i,dt,dphimax)
!       'i'-loop should be inner loop because phi(i,k) is contiguous with index i
!$OMP   DO REDUCTION(max:dphimax)
        do k=1,kmax-1
          do i=1,imax-1
            dphi=(phi(i+1,k)+phi(i-1,k)-2.*phi(i,k))*dy2i
     f           +(phi(i,k+1)+phi(i,k-1)-2.*phi(i,k))*dx2i
            dphi=dphi*dt
            dphimax=max(dphimax,dphi)
            phin(i,k)=phi(i,k)+dphi
          enddo
        enddo
!$OMP   END DO
! save values
!       'i'-loop should be inner loop because phi(i,k) is contiguous with index i
!$OMP   DO
        do k=1,kmax-1
          do i=1,imax-1
            phi(i,k)=phin(i,k)
          enddo
        enddo
!$OMP   END DO
!$OMP  END PARALLEL
       if(dphimax.lt.eps) goto 10
      enddo
10    continue
!
      call cpu_time(t1)
!$    wt1=omp_get_wtime() 
! print array
      write (*,'(/,a,i6,a)')
     f      'phi after',it,' iterations'
      call heatpr(phi,imax,kmax)
      write (*,'(/,a,1pg12.4)')  'cpu_time  :   ', t1-t0
!$    write (*,'(/,a,1pg12.4)')  'omp_get_wtime:', wt1-wt0
!
      stop
      end
!
!
!
!
      subroutine heatpr(phi,imx,kmx)
      double precision phi(0:imx,0:kmx)
      integer imx, kmx 
!
      integer i, prt_i,prt_k, prt_i_max, prt_k_max
      integer k(20) 

      prt_i_max = min(prtrows, imx+1) 
      prt_k_max = min(prtcols, kmx+1, 20) 
      do prt_k=0, prt_k_max-1
        k(prt_k) = 1.0*prt_k/(prt_k_max-1)*kmx + 0.5/*rounding*/
      end do
      write (*,'(''i  |  k='',20(''  '',I6))')
     +         (k(prt_k), prt_k=0,prt_k_max-1)
      write (*,'(''-------+'',20A)') ('--------', prt_k=0,prt_k_max-1)
      do prt_i=0, prt_i_max-1
        i=1.0*prt_i/(prt_i_max-1)*imx + 0.5/*rounding*/
        write (*, '(I6,'' |'',20(2X,F6.4))') 
     +              i, (phi(i,k(prt_k)), prt_k=0,prt_k_max-1) 
      end do
!
      return
      end
