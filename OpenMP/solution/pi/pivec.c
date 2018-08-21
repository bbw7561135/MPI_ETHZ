/* pivec.c - This is a demonstration program to show
             cache false-sharing between OpenMP threads */

/* CAUTION: This program is not intended to be the
            basis of an efficient OpenMP parallelization. */

/* Compilation:
     without -DTOTAL :
       The data is written with static schedule with a first 
       loop into a vector array and 2 times, half of the vector
       is fetched in a second and a third loop.
     with -DTOTAL :
       The data is written and read with the same vector size.
       To achieve similar loop overhead, a short third loop is added.
   Execution:
     The first argument defines the vector size (in doubles).
*/ 

/* Expected result: As long as m is less than the cache size,
   there should be a significant panalty if -DTOTAL is not used.
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
#include <stdlib.h>

#define f(A) (4.0/(1.0+A*A))

const int n = 10000000;
      int m = 10000000; /*vector chunk*/

int main(int argc, char** argv)
{
  int i,j,j_e;
  double w,x,sum,pi,sum0, *vec;
  clock_t t1,t2;
  struct timeval tv1,tv2; struct timezone tz;
# ifdef _OPENMP
    double wt1,wt2;
# endif
 
  if (argc>1) m=atoi(argv[1]);
  vec = (double *)malloc(m*sizeof(double));
 
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
for (j=1;j<=n;j=j+m)
{
 j_e = (j+m-1>n ? n : j+m-1);
#pragma omp parallel private (x,sum0) shared(w,sum)
 {
  sum0=0;
#pragma omp for schedule(static)
  for (i=j;i<=j_e;i++)
  {
    x=w*((double)i-0.5);
    vec[i-j] = f(x);
  }

#ifdef TOTAL

#pragma omp for schedule(static)
  for (i=j;i<=j_e;i++)
  {
    sum0=sum0+vec[i-j];
  }
#pragma omp for schedule(static)
  for (i=0;i<16;i++)
  {
    sum0=1.0*sum0;
  }

#else

#pragma omp for schedule(static)
  for (i=j;i<(j+j_e)/2;i++)
  {
    sum0=sum0+vec[i-j];
  }
#pragma omp for schedule(static)
  for (i=(j+j_e)/2;i<=j_e;i++)
  {
    sum0=sum0+vec[i-j];
  }

#endif

#pragma omp critical
  {
    sum=sum+sum0;
  }
 }
} /*end loop j*/
  pi=w*sum;

  t2=clock();
# ifdef _OPENMP
    wt2=omp_get_wtime();
# endif
  gettimeofday(&tv2, &tz);
  printf( "computed pi = %24.16g\n", pi );
  printf( "vector chunk size = %d doubles = %d bytes\n", m, m*sizeof(double));
  printf( "CPU time (clock)                = %12.4g sec\n", (t2-t1)/1000000.0 );
# ifdef _OPENMP
    printf( "wall clock time (omp_get_wtime) = %12.4g sec\n", wt2-wt1 );
# endif
  printf( "wall clock time (gettimeofday)  = %12.4g sec\n", (tv2.tv_sec-tv1.tv_sec) + (tv2.tv_usec-tv1.tv_usec)*1e-6 );

  return 0;
}
