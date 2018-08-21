#include <stdio.h>

# ifdef _OPENMP
#  include <omp.h>
# endif

int main(int argc, char** argv)
{
  int i;

  i = -1;

# ifdef _OPENMP
  i = omp_get_thread_num();
# endif

  printf( "hello world %d\n", i );
  return 0;
}
