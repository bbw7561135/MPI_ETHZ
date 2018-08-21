program compute_pi
implicit none 

! times using cpu_time
real t0,t1
!--unused-- use omp_lib 
!$ integer omp_get_num_threads 
!$ double precision omp_get_wtime
!$ double precision wt0,wt1   

integer i
integer, parameter        :: n=10000000
real(kind=8) w,x,sum,pi,f,a
integer i_outer 
real(kind=8) sum_outer(0:63,0:99) ! A cache line should be not larger than the 64 ints.
                                  ! Each task j should use sum_outer(0,j), i.e., each of
                                  ! these summing-variables is in a different cache-line

! function to integrate
f(a)=4.0_8/(1.0_8+a*a)

!$omp parallel
!$omp single
!$  write (*,*) 'OpenMP-parallel with',omp_get_num_threads(),'threads'
!$omp end single
!$omp end parallel

!$ wt0=omp_get_wtime() 
call cpu_time(t0)

! calculate pi = integral [0..1] 4/(1+x**2) dx
w=1.0_8/n
sum=0.0_8

do i_outer=0,99 
  sum_outer(0,i_outer)=0.0
enddo

!$OMP PARALLEL
!$OMP  SINGLE
!$OMP   TASKLOOP PRIVATE(x,i)
  do i_outer=0,99
    do i=i_outer+1,n,100
      x=w*(i-0.5_8)
      sum_outer(0,i_outer) =  sum_outer(0,i_outer) + f(x);
    enddo
  enddo
!$OMP   END TASKLOOP
!$OMP  END SINGLE
!$OMP END PARALLEL

do i_outer=0,99 
  sum = sum + sum_outer(0,i_outer)
enddo
pi=w*sum

call cpu_time(t1)
!$ wt1=omp_get_wtime() 

write (*,'(/,a,1pg24.16)')   'computed pi = ', pi
write (*,'(/,a,1pg12.4)')    'cpu_time:     ', t1-t0
!$ write (*,'(/,a,1pg12.4)') 'omp_get_wtime:', wt1-wt0
stop
end program compute_pi
