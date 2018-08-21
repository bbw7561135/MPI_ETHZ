program confl

integer, parameter:: N=10
real:: a(0:N), b(N)

a(0) = 0.

!$OMP PARALLEL

!$OMP DO
DO i=1,N
  a(i) = 1./i
  b(i) = a(i-1)+a(i)
ENDDO

!$OMP END PARALLEL

print*, b

end program confl
