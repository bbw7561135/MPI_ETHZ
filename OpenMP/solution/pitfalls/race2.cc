#include <iostream>

int main(void){
  int jlast=-777;
  int i;
#pragma omp parallel for lastprivate(jlast)
  for(i=1; i<1000; i++){
    jlast=i;
  }

  std::cout << jlast << "\n";


  return 0;
}
