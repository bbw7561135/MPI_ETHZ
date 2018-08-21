#
set xlabel 'size of loop'
set ylabel 'Mops/s'
set log
set nolog y



set title 'OpenMP norm benchmark C'
plot \
'solution_c.log' u 1:4 t 'solution' w lp, \
'result.log' u 1:4 t 'my result' w lp
pause -1


