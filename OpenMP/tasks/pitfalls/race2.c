#include <stdio.h>

int main(void){
  int jlast=-777;
  int i;
#pragma omp parallel for private(jlast)
  for(i=1; i<1000; i++){
    jlast=i;
  }

  printf("%d\n",jlast);


  return 0;
}
