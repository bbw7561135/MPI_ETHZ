
real function WORK1(A)
integer A

WORK1=real(A+1)
end function WORK1


real function WORK2(A,B)
real A
integer B

WORK2=A+real(B)
end function WORK2

program race1
implicit none


real WORK1, WORK2
real X, TMP
real Y
integer I

       X=0
!$OMP PARALLEL SHARED (X), PRIVATE(TMP) 

!$OMP DO REDUCTION(+:X)
       DO 100 I=1,100
              TMP = WORK1(I)
              X = X + TMP
100  CONTINUE
!$OMP END DO NOWAIT

       Y = X
!$OMP END PARALLEL   

print *,Y



end program race1
