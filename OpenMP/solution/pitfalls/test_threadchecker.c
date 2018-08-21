#include <stdio.h>
#ifdef _OPENMP
#  include <omp.h>
#endif
/* for sleep(): */
#include <unistd.h>

const int n = 4;

int main(int argc, char** argv)
{
  int i, sum, psum;
  int thread_num[n], num_threads[n], max_threads[n], num_procs[n], in_parallel[n], dynamic[n], nested[n];
  double wtime[n], wtick[n];
  int a[n];

# ifndef _OPENMP
    printf("This program is not compiled with OpenMP\n");
# endif

#   pragma omp parallel
    {
#     pragma omp for nowait schedule(static,2)
        for (i=0; i<n; i++) 
		{
		  a[i] = 2*i;
		}
#     pragma omp for schedule(static,2)
        for (i=0; i<n; i++) 
		{
		  a[i] = a[i] + 3;
		}
#     pragma omp for schedule(static,2)
        for (i=0; i<n; i++) 
		{
		  printf("a[%1d]=%2d  ",i,a[i]);
		}
		printf("\n");
    } /* end omp parallel */
      
 
    { 
      int thread_num, num_threads, max_threads, num_procs, in_parallel, dynamic, nested;
      double wtime, wtick;
      thread_num  = omp_get_thread_num(); 
      num_threads = omp_get_num_threads();
      max_threads = omp_get_max_threads();
      num_procs   = omp_get_num_procs();
      in_parallel = omp_in_parallel();
      dynamic     = omp_get_dynamic();
      nested      = omp_get_nested(); 
	  sleep(2);
      wtime       = omp_get_wtime(); 
      wtick       = omp_get_wtick(); 
      printf("I am thread_num %2d of %2d num_threads, %2d max_threads, %2d num_procs, %2d in_parallel, %2d dynamic, %2d nested, %e wtime, %e wtick - OUTSIDE PARALLEL REGION\n", 
	               thread_num,        num_threads,  max_threads,  num_procs,  in_parallel,  dynamic,  nested, wtime, wtick); 
    }

	omp_set_num_threads(4);
	omp_set_dynamic(1);
	omp_set_nested(1);
	
    { 
      int thread_num, num_threads, max_threads, num_procs, in_parallel, dynamic, nested;
      double wtime, wtick;
      thread_num  = omp_get_thread_num(); 
      num_threads = omp_get_num_threads();
      max_threads = omp_get_max_threads();
      num_procs   = omp_get_num_procs();
      in_parallel = omp_in_parallel();
      dynamic     = omp_get_dynamic();
      nested      = omp_get_nested(); 
	  sleep(2);
      wtime       = omp_get_wtime(); 
      wtick       = omp_get_wtick(); 
      printf("I am thread_num %2d of %2d num_threads, %2d max_threads, %2d num_procs, %2d in_parallel, %2d dynamic, %2d nested, %e wtime, %e wtick - OUTSIDE PARALLEL REGION\n", 
	               thread_num,        num_threads,  max_threads,  num_procs,  in_parallel,  dynamic,  nested, wtime, wtick); 
    }

    sum = 0;
#   pragma omp parallel 
    { 
      int thread_num, num_threads, max_threads, num_procs, in_parallel, dynamic, nested;
      double wtime, wtick;
      thread_num  = omp_get_thread_num(); 
      num_threads = omp_get_num_threads();
      max_threads = omp_get_max_threads();
      num_procs   = omp_get_num_procs();
      in_parallel = omp_in_parallel();
      dynamic     = omp_get_dynamic();
      nested      = omp_get_nested(); 
	  sleep(2);
      wtime       = omp_get_wtime(); 
      wtick       = omp_get_wtick(); 
      printf("I am thread_num %2d of %2d num_threads, %2d max_threads, %2d num_procs, %2d in_parallel, %2d dynamic, %2d nested, %e wtime, %e wtick - INSIDE PARALLEL REGION\n", 
	               thread_num,        num_threads,  max_threads,  num_procs,  in_parallel,  dynamic,  nested, wtime, wtick); 
#     pragma omp critical
      {
       printf("I am thread_num %2d of %2d num_threads, %2d max_threads, %2d num_procs, %2d in_parallel, %2d dynamic, %2d nested, %e wtime, %e wtick - INSIDE CRITICAL\n", 
	               thread_num,        num_threads,  max_threads,  num_procs,  in_parallel,  dynamic,  nested, wtime, wtick); 
       sum = sum + 7;
      } /* end omp critical */   
    } /* end omp parallel */
	printf("sum = %2d \n", sum);

    sum = 0;
#   pragma omp parallel private(psum)
    { 
#     pragma omp for
        for (i=0; i<n; i++)
		{
		  psum = 7;
		}  /* end omp for */
#     pragma omp critical
      {
       sum = sum + psum;
      } /* end omp critical */   
    } /* end omp parallel */
	printf("sum = %2d \n", sum);

    sum = 0;
#   pragma omp parallel private(psum)
    { 
#     pragma omp for schedule(static,1)
        for (i=0; i<n; i++)
		{
		  psum = 7;
		}  /* end omp for */
#     pragma omp critical
      {
       sum = sum + psum;
      } /* end omp critical */   
    } /* end omp parallel */
	printf("sum = %2d \n", sum);

    for (i=0; i<n; i++)
    { 
      thread_num[i]  = -111;
      num_threads[i] = -999;
	  max_threads[i] = -999;
	  num_procs[i]   = -999;
	  in_parallel[i] = -999;
	  dynamic[i]     = -999;
	  nested[i]      = -999;
	  wtime[i]       = -999;
	  wtick[i]       = -999;
    } 
 
#   pragma omp parallel for 
    for (i=0; i<n; i++)
    { 
      thread_num[i]  = /*omp_get_thread_num()*/ 0; 
      num_threads[i] = omp_get_num_threads();
      max_threads[i] = omp_get_max_threads();
      num_procs[i]   = omp_get_num_procs();
      in_parallel[i] = omp_in_parallel();
      dynamic[i]     = omp_get_dynamic();
      nested[i]      = omp_get_nested(); 
	  sleep(i+1);
      wtime[i]       = omp_get_wtime(); 
      wtick[i]       = omp_get_wtick(); 
      printf("I am thread_num %2d of %2d num_threads, %2d max_threads, %2d num_procs, %2d in_parallel, %2d dynamic, %2d nested, %e wtime, %e wtick - INSIDE LOOP\n", 
	               thread_num[i],        num_threads[i],  max_threads[i],  num_procs[i],  in_parallel[i],  dynamic[i],  nested[i], wtime[i], wtick[i]); 
    } /* end omp parallel */

    for (i=0; i<n; i++)
    { 
      printf("I am thread_num %2d of %2d num_threads, %2d max_threads, %2d num_procs, %2d in_parallel, %2d dynamic, %2d nested, %e wtime, %e wtick \n", 
	               thread_num[i],        num_threads[i],  max_threads[i],  num_procs[i],  in_parallel[i],  dynamic[i],  nested[i], wtime[i], wtick[i]); 
    } /* end omp parallel */

  return 0;
}
