      program compute_pi
      implicit none 
      integer n,i
      double precision w,x,sum,pi,f,a
      parameter (n=10 000 000)
! times using cpu_time
      real t0
      real t1
!--unused-- include 'omp_lib.h' 
!$    double precision omp_get_wtime
!$    double precision wt0,wt1   
!
!$    integer omp_get_thread_num, omp_get_num_threads 
!$    integer myrank, num_threads
      logical openmp_is_used
!
! function to integrate
      f(a)=4.d0/(1.d0+a*a)
!
!$omp parallel private(myrank, num_threads)
!$      myrank = omp_get_thread_num()
!$      num_threads = omp_get_num_threads()
!$      write (*,*) 'I am thread',myrank,'of',num_threads,'threads' 
!$omp end parallel
      openmp_is_used = .false. 
!$    openmp_is_used = .true. 
      if (.not. openmp_is_used) then
        write (*,*) 'This program is not compiled with OpenMP'
      endif
!
!$    wt0=omp_get_wtime() 
      call cpu_time(t0)
! 
! calculate pi = integral [0..1] 4/(1+x**2) dx
      w=1.0d0/n
      sum=0.0d0
      do i=1,n
        x=w*(i-0.5d0)
        sum=sum+f(x)
      enddo
      pi=w*sum
!
      call cpu_time(t1)
!$    wt1=omp_get_wtime() 
      write (*,'(/,a,1pg24.16)') 'computed pi = ', pi
      write (*,'(/,a,1pg12.4)')  'cpu_time  :   ', t1-t0
!$    write (*,'(/,a,1pg12.4)')  'omp_get_wtime:', wt1-wt0
      stop
      end
