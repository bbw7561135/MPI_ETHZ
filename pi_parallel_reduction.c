#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#ifdef _OPENMP
#  include <omp.h>
#endif

#define f(A) (4.0/(1.0+A*A))

const int n = 10000000;

int main(int argc, char** argv)
{
  int i;
  double w,x,sum,pi;
  clock_t t1,t2;
  struct timeval tv1,tv2; struct timezone tz;
# ifdef _OPENMP
    double wt1,wt2;
# endif

/* Begin of SPACE for FIRST EXERCISE */
# ifdef _OPENMP
	int total_thread;
	int thread_num;
	/*private(thread_num)*/
	/*# pragma omp parallel private(thread_num)*/
	# pragma omp parallel private(thread_num, total_thread)
	{
		total_thread = omp_get_num_threads();
		thread_num = omp_get_thread_num();
		/*sleep(1);*/
		sleep(1);
		printf("I am thread %i of %i threads.\n", thread_num, total_thread);
	}
# else
	printf("This program is not complied with OpenMP.\n");
# endif


/* End of SPACE for FIRST EXERCISE */

  gettimeofday(&tv1, &tz);
# ifdef _OPENMP
    wt1=omp_get_wtime();
# endif
  t1=clock();
 
/* calculate pi = integral [0..1] 4/(1+x**2) dx */
w=1.0/n;
sum=0.0; 
double psum; 
# ifdef _OPENMP
	//# pragma omp parallel private(w, sum)
	//{
	//w=1.0/n;
	//sum=0.0;
	#pragma omp parallel for private(x), shared(w), reduction(+:sum)
	for (i=1;i<=n;i++)
	{
		x=w*((double)i-0.5);
		//# pragma omp critical
		sum=sum+f(x);
	}
	# pragma omp critical
	//sum += psum;
	/*# pragma omp parallel private(psum), shared(w, sum)
	{
		psum = 0.0;
		# pragma omp for private (x)
		for (i=1;i<=n;i++)
		{
			x=w*((double)i-0.5);
			//# pragma omp critical
			psum=psum+f(x);
		}
		# pragma omp critical
		sum += psum;
	}*/
	/*# pragma omp for private (x)
	for (i=1;i<=n;i++)
	{
		x=w*((double)i-0.5);
		//# pragma omp critical
		sum=sum+f(x);
	}
	# pragma omp critical
	//printf("x is %f.\n", x);*/
	pi=w*sum;
	//}
# else
	w=1.0/n;
	sum=0.0;
	for (i=1;i<=n;i++)
	{
		x=w*((double)i-0.5);
		sum=sum+f(x);
	}
	pi=w*sum;
# endif
 
  t2=clock();
# ifdef _OPENMP
    wt2=omp_get_wtime();
# endif
  gettimeofday(&tv2, &tz);
  printf( "computed pi = %24.16g\n", pi );
  printf( "CPU time (clock)                = %12.4g sec\n", (t2-t1)/1000000.0 );
# ifdef _OPENMP
    printf( "wall clock time (omp_get_wtime) = %12.4g sec\n", wt2-wt1 );
# endif
  printf( "wall clock time (gettimeofday)  = %12.4g sec\n", (tv2.tv_sec-tv1.tv_sec) + (tv2.tv_usec-tv1.tv_usec)*1e-6 );
  return 0;
}
