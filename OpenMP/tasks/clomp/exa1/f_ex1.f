!This simple example shows that for some situations, 
!no changes are necessary to an existing OpenMP code.
      include 'omp_lib.h'
      
!These two variables are autopromoted to sharable
!by the compiler
      integer a, num_threads
      a = 0

!$omp parallel 
!$omp single 
	  num_threads = omp_get_num_threads()
!$omp end single nowait

!$omp critical
      a = a + 1;
!$omp end critical

!$omp end parallel 

20     format ("a=", I3, ".  It should be ", I3, ".")
       write (6,20), a, num_threads
       end
