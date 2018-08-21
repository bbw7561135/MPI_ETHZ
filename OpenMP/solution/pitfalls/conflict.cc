
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

const int NMAX=11;

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
  }

#pragma omp for
  for(i=1; i<NMAX; i++){
    b[i]=a[i-1]+a[i];
  }
  }/* end parallel*/

  for(i=1;i<NMAX; i++){
    std::cout << b[i] << "\n";
  }

  return 0;
}
