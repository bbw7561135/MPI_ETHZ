
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define NMAX 11

double a[NMAX];
double b[NMAX];


int main(void){
  int i;
  a[0]=0.0;

#pragma omp parallel 
  {
#pragma omp for
  for(i=1; i<NMAX; i++){
    a[i]=1.0/i;
	if (i == 1){
		b[i]=a[i];
	}
	else
	{
		b[i]=1.0/(i -1)+a[i];
	}
    
  }
  }/* end parallel*/

  for(i=1;i<NMAX; i++){
    printf("%g\n",b[i]);
  }

  return 0;
}
