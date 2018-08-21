set term post color
#
set xlabel 'size of loop'
set ylabel 'Mops/s'
set log
set nolog y
set out 'result-norm.eps'


set title 'OpenMP norm benchmark Fortran'
plot 'result-norm_f90.log' u 1:4 t 'single' w lp, \
'result-norm_f90.log' u 1:(2*$4) t 'speed up 2' w l, \
'result-norm_f90_omp.log' u 1:4 t 'OpenMP' w lp 

set title 'OpenMP norm benchmark C'
plot 'result-norm_c.log' u 1:4 t 'single' w lp, \
'result-norm_c.log' u 1:(2*$4) t 'speed up 2' w l, \
'result-norm_c_omp.log' u 1:4 t 'OpenMP' w lp 

set title 'OpenMP norm benchmark PC'
plot 'pc/result-norm_f90.log' u 1:4 t 'single' w lp, \
'pc/result-norm_f90.log' u 1:(2*$4) t 'speed up 2' w l, \
'pc/result-norm_f90_omp.log' u 1:4 t 'OpenMP' w lp ,\
'pc/result-norm_f90_omp_if.log' u 1:4 t 'OpenMP IF' w lp

set title 'OpenMP norm benchmark SX-5'
plot 'sx5/result-norm_f90.log' u 1:4 t 'single' w lp, \
'sx5/result-norm_f90.log' u 1:(2*$4) t 'speed up 2' w l, \
'sx5/result-norm_f90_omp.log' u 1:4 t 'OpenMP' w lp

set title 'OpenMP norm benchmark HPV'
plot 'hpv/result-norm_f90.log' u 1:4 t 'single' w lp, \
'hpv/result-norm_f90_1thread.log' u 1:4 t '1 thread' w l, \
'hpv/result-norm_f90.log' u 1:(2*$4) t 'speed up 2' w l, \
'hpv/result-norm_f90_omp.log' u 1:4 t 'OpenMP' w lp

set title 'OpenMP norm benchmark SGI'
plot 'sgi/result-norm_f90.log' u 1:4 t 'single' w lp, \
'sgi/result-norm_f90.log' u 1:(2*$4) t 'speed up 2' w l, \
'sgi/result-norm_f90_omp.log' u 1:4 t 'OpenMP' w lp

set title 'OpenMP norm benchmark SR8K'
plot 'sr8k/result-norm_f90.log' u 1:4 t 'single' w lp, \
'sr8k/result-norm_f90.log' u 1:(2*$4) t 'speed up 2' w l, \
'sr8k/result-norm_f90_omp.log' u 1:4 t 'OpenMP' w lp

set title 'OpenMP vs. Compas SR8K'
plot 'sr8k/result-norm_f90.log' u 1:4 t 'single' w lp, \
'sr8k/result-norm_f90.log' u 1:(8*$4) t 'speed up 8' w l, \
'sr8k/result-norm_f90_omp_8.log' u 1:4 t 'OpenMP' w lp,\
'sr8k/result-norm_f90_compass.log' u 1:4 t 'compass' w lp

set title 'scalar norm benchmark'
plot 'sr8k/result-norm_f90.log' u 1:4 t 'sr8k' w lp, \
'sr8k/result-norm_f90_compass.log' u 1:4 t 'sr8k compass' w lp, \
 'pc/result-norm_f90.log' u 1:4 t 'pc' w lp, \
 'sgi/result-norm_f90.log' u 1:4 t 'sgi' w lp, \
 'hpv/result-norm_f90.log' u 1:4 t 'hpv' w lp, \
 'sx5/result-norm_f90.log' u 1:4 t 'sx5' w lp

set title 'OpenMP norm benchmark'
plot 'sr8k/result-norm_f90_omp.log' u 1:4 t 'sr8k' w lp, \
 'pc/result-norm_f90_omp.log' u 1:4 t 'pc' w lp, \
 'sgi/result-norm_f90_omp.log' u 1:4 t 'sgi' w lp, \
 'hpv/result-norm_f90_omp.log' u 1:4 t 'hpv' w lp, \
 'sx5/result-norm_f90_omp.log' u 1:4 t 'sx5' w lp




