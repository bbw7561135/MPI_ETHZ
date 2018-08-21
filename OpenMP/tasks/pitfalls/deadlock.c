#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define N 10


int main(void){
  double x,y;
  int i;
  omp_lock_t lock_x, lock_y;
  
  omp_init_lock(&lock_x);
  omp_init_lock(&lock_y);
  x = 0;
  y = 0;

# pragma omp parallel shared(x,y)
  {
#   pragma omp for 
    for(i=1; i<=N; i++){
	  if (i < 0.3*N) {
	    omp_set_lock(&lock_x);
	      x = x + i;
	      omp_set_lock(&lock_y);
	        y = y + i;
	      omp_unset_lock(&lock_y);
	    omp_unset_lock(&lock_x);
	  } else {
	    omp_set_lock(&lock_y);
	      y = y + i;
	      omp_set_lock(&lock_x);
	        x = x + i;
	      omp_unset_lock(&lock_x);
	    omp_unset_lock(&lock_y);
	  } /*endif*/
    }

  }/* end parallel*/
  printf("x=%g, y=%g, expected=%g, on %i threads\n",x,y,(N+1.0)*N/2.0,omp_get_max_threads());
}
