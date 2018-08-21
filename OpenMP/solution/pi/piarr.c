/* piarr.c - This is a demonstration program to show
             cache-line false-sharing between OpenMP threads */

/* CAUTION: This program is not intended to be the
            basis of an efficient OpenMP parallelization. */

/* Execution:
     The first argument defines the stride (in doubles).
     The value must be between 1 and MAX_STRIDE (defined below). 
     If the program should be executed with more threads than
     defined in MAX_THREADS, MAX_THREADS must be modified.
*/ 

/* Expected result: 
     Each thread accesses always the same array-element.
     If the stride is larger than the cache-line size
     then the performance should be similar to that of pic2.c,
     except that always memory accesses are done.
     If the stride is less than the the cache-line size
     then a very slow execution is expected.  
*/
/* Author: Rolf Rabenseifner, HLRS, University of Stuttgart, June 2005
           http://www.hlrs.de/people/rabenseifner/  
*/
 
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#ifdef _OPENMP
#  include <omp.h>
#endif

#define MAX_THREADS 16
#define MAX_STRIDE  100 /* in doubles */

#define f(A) (4.0/(1.0+A*A))

const int n = 100000000;

int main(int argc, char** argv)
{
  int i,myrank,stride;
  double w,x,sum,pi,p_sum[MAX_THREADS*MAX_STRIDE];
  clock_t t1,t2;
  struct timeval tv1,tv2; struct timezone tz;
# ifdef _OPENMP
    double wt1,wt2;
# endif

  stride=1;
  if (argc>1) stride=atoi(argv[1]);

# ifdef _OPENMP
#   pragma omp parallel
    { 
#     pragma omp single 
      printf("OpenMP-parallel with %1d threads\n", omp_get_num_threads());
    } /* end omp parallel */
# endif

  gettimeofday(&tv1, &tz);
# ifdef _OPENMP
    wt1=omp_get_wtime();
# endif
  t1=clock();
  
/* calculate pi = integral [0..1] 4/(1+x**2) dx */
  w=1.0/n;
  sum=0.0;
#pragma omp parallel private (x,myrank) shared(w,sum,p_sum)
{
# ifdef _OPENMP
    myrank=stride*omp_get_thread_num();
# else
    myrank=0;
# endif
  p_sum[myrank]=0;
#pragma omp for
  for (i=1;i<=n;i++)
  {
    x=w*((double)i-0.5);
    p_sum[myrank+(x>1?1:0)] = p_sum[myrank] + f(x);
    /* The term (x>1?1:0) is always zero. It is used to prohibit 
       register caching of of p_sum[myrank], i.e., to guarantee that
       each access to this variable is done via cache in the memory. */
  }
#pragma omp critical
  {
    sum=sum+p_sum[myrank];
  }
}
  pi=w*sum;

  t2=clock();
# ifdef _OPENMP
    wt2=omp_get_wtime();
# endif
  gettimeofday(&tv2, &tz);
  printf( "computed pi = %24.16g\n", pi );
  printf( "stride      = %24d\n",    stride );
  printf( "CPU time (clock)                = %12.4g sec\n", (t2-t1)/1000000.0 );
# ifdef _OPENMP
    printf( "wall clock time (omp_get_wtime) = %12.4g sec\n", wt2-wt1 );
# endif
  printf( "wall clock time (gettimeofday)  = %12.4g sec\n", (tv2.tv_sec-tv1.tv_sec) + (tv2.tv_usec-tv1.tv_usec)*1e-6 );
  return 0;
}
