!
!	Computing matrix norm in parallel using OpenMP
!	Prints 
!		size of matrix
!		real time
!		user time
!		number of iterations (for timimg loop)
!
!	V2 using, only one parallel section
!
!	by Holger Berger/NEC ESS 
!       and Matthias Mueller/HLRS University Stuttgart
!
program norm_test

	double precision:: totals,reals,users
	double precision:: totale,reale,usere
	real :: s

	real:: a(4097,4097)

	integer :: n,i,steps
	integer :: cpus, threads

!$	integer OMP_GET_NUM_THREADS
!$	integer OMP_GET_NUM_PROCS
        cpus=0
        threads=0
!$OMP PARALLEL
!$OMP MASTER
!$ 	cpus = OMP_GET_NUM_PROCS()
!$ 	threads = OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL

        call hlrs_get_time(usere,totals)


	print *, '# OpenMP Version, ',threads,' threads on ',cpus,' CPUs '
	print *, '# size   real  repetitions Mops/s'

	a = 1.0

        print *, '# Triangular Matrix: '
	n=2
	do steps=1,12
		i = 0
		call norm(s,a,n)
		call hlrs_get_time(users,reals)
		do
			i = i+1
			call norm(s,a,n)
			call hlrs_get_time(usere,reale)
			if (usere-users > 0.2 .or. i>100000000) exit
		end do
		print *, n,(reale-reals)/dble(i),i,dble(i)*n*n/(reale-reals)
		if (s .ne. n) print *,'# wrong result!'
		n=n*2
	end do

	n=2
	do steps=1,12
		i = 0
		call norm_triangular(s,a,n)
		call hlrs_get_time(users,reals)
		do
			i = i+1
			call norm_triangular(s,a,n)
			call hlrs_get_time(usere,reale)
			if (usere-users > 0.2 .or. i>100000000) exit
		end do
		print *,"#", n,(reale-reals)/dble(i),i,dble(i)*n*n/(reale-reals)
		if (s .ne. n) print *,'# wrong result!'
		n=n*2
	end do
        call hlrs_get_time(usere,totale)
        print *,"# total time = ",totale-totals, " seconds"
end program norm_test

subroutine norm(s,a,n)
	real :: s,a(n,n),mymax
	integer :: n,i,j
        s=0

!$OMP PARALLEL DO PRIVATE(i,mymax) REDUCTION(MAX:s) 
	do j=1,n
	  mymax=0
	  do i=1,n
	    mymax = mymax + abs(a(i,j))
	  end do
          s=max(s,mymax)
 	end do




end subroutine norm



subroutine norm_triangular(s,a,n)
	real :: s,a(n,n),mymax
	integer :: n,i,j
        s=0

!$OMP PARALLEL DO PRIVATE(i,mymax) REDUCTION(MAX:s) SCHEDULE(GUIDED)
	do j=1,n
	  mymax=0
	  do i=1,j
	    mymax = mymax + abs(a(i,j))
	  end do
          s=max(s,mymax)
 	end do




end subroutine norm_triangular
