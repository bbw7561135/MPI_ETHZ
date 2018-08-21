/* This simple example shows that for some situations, 
 * no changes are necessary to an existing OpenMP code.
 */
#include <omp.h>

int main()
{
   /* These two variables are autopromoted to sharable 
    * by the compiler */
   int a = 0; 
   int num_threads; 

#pragma omp parallel 
   {
      #pragma omp single nowait
	  num_threads = omp_get_num_threads();

      #pragma omp critical
      a++;
   }

   printf("a=%i.  It should be %i.\n", a, num_threads);
}
