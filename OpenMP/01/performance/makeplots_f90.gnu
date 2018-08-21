#
set xlabel 'size of loop'
set ylabel 'Mops/s'
set log
set nolog y



set title 'OpenMP norm benchmark Fortran'
plot \
'solution_f90.log' u 1:4 t 'solution' w lp, \
'result.log' u 1:4 t 'my result' w lp
pause -1


