program calc1
      JLAST = -777
!$OMP PARALLEL DO PRIVATE(JLAST)
      DO 10 J=1,1000
           JLAST = J
10    CONTINUE 
!$OMP END PARALLEL DO
      print *, JLAST
end program calc1
