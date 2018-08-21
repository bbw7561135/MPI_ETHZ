program compute_pi
implicit none 

! times using cpu_time
real t0,t1
!--unused-- use omp_lib 
!$ double precision omp_get_wtime
!$ double precision wt0,wt1   

integer i
integer, parameter        :: n=10000000
real(kind=8) w,x,sum,pi,f,a

! function to integrate
f(a)=4.0_8/(1.0_8+a*a)

! Begin of SPACE for FIRST EXERCISE


! End of SPACE for FIRST EXERCISE

!$ wt0=omp_get_wtime() 
call cpu_time(t0)

! calculate pi = integral [0..1] 4/(1+x**2) dx
w=1.0_8/n
sum=0.0_8
do i=1,n
   x=w*(i-0.5_8)
   sum=sum+f(x)
enddo
pi=w*sum

call cpu_time(t1)
!$ wt1=omp_get_wtime() 

write (*,'(/,a,1pg24.16)')   'computed pi = ', pi
write (*,'(/,a,1pg12.4)')    'cpu_time:     ', t1-t0
!$ write (*,'(/,a,1pg12.4)') 'omp_get_wtime:', wt1-wt0
stop
end program compute_pi
