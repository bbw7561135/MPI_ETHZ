! CAUTION: This program is not intended to be the
!          basis of an efficient OpenMP parallelization. */
!
      program compute_pi
      implicit none 
      integer n,i
      double precision w,x,sum,myf,pi,f,a
      parameter (n=10 000 000)
! times using cpu_time
      real t0
      real t1
!--unused-- include 'omp_lib.h' 
!$    integer omp_get_num_threads 
!$    double precision omp_get_wtime
!$    double precision wt0,wt1   
!
! function to integrate
      f(a)=4.d0/(1.d0+a*a)
!
!$omp parallel
!$omp single
!$    write(*,*)'OpenMP-parallel with',omp_get_num_threads(),'threads'
!$omp end single
!$omp end parallel
!
!$    wt0=omp_get_wtime() 
      call cpu_time(t0)
! 
! calculate pi = integral [0..1] 4/(1+x**2) dx
      w=1.0d0/n
      sum=0.0d0
!$OMP PARALLEL PRIVATE(x,myf), SHARED(w,sum)
!$OMP DO ORDERED
      do i=1,n
        x=w*(i-0.5d0)
        myf=f(x) 
!       ! The ordered region inside of a loop with many
!       ! iteration may cause a huge overhead 
!$OMP   ORDERED
          sum=sum+myf
!$OMP   END ORDERED
      enddo
!$OMP END DO
!$OMP END PARALLEL
      pi=w*sum
!
      call cpu_time(t1)
!$    wt1=omp_get_wtime() 
      write (*,'(/,a,1pg24.16)') 'computed pi = ', pi
      write (*,'(/,a,1pg12.4)')  'cpu_time  :   ', t1-t0
!$    write (*,'(/,a,1pg12.4)')  'omp_get_wtime:', wt1-wt0
      stop
      end
