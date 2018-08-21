      module dims
      real*8,allocatable :: xdim(:)
      real*8,allocatable :: ydim(:)
      end module

      use dims
      include 'omp_lib.h'
      integer n, niter
      data n/2000/, niter/10/
      integer i
      real*8 pot

      allocate(xdim(n))
      allocate(ydim(n))

      call initPositions(n)
      call updatePositions(n)

      do i=1, niter
         pot = computePot(n)
	 print *, "Potential: ", pot
	 call updatePositions(n)
      enddo

      end

      subroutine initPositions(n)
      use dims
      integer i, j
      integer n
      real*8 num

      do j=1,n
	 call random_number(num)
	 xdim(j) = 0.5 + num
	 call random_number(num)
	 ydim(j) = 0.5 + num
      enddo

      end subroutine

      subroutine updatePositions(n)
      use dims
      integer i, j
      integer n
      real*8 num

      do j=1,n
	 call random_number(num)
	 xdim(j) = xdim(j) - (0.5 + num)
	 call random_number(num)
	 ydim(j) = ydim(j) - (0.5 + num)
      enddo

      end subroutine


      real*8 function computePot(n) result(pot)
      integer i, j
      integer n
      real*8 distx,disty,dist
       pot = 0

!$omp parallel private(distx, disty, dist, i, j) 
!$omp do schedule(static) reduction(+:pot)
      do i=1,n
	  do j=1,i-1
	     call power(distx, disty, i, j)
	     dist = sqrt(distx + disty)
	     pot = pot + (1.0 / dist)
	  enddo
      enddo
!$omp end parallel

      end function

      subroutine power(distx, disty, i, j)
      use dims
      real*8 distx, disty
      integer i, j

      distx = (xdim(j) - xdim(i)) * (xdim(j) - xdim(i))
      disty = (ydim(j) - ydim(i)) * (ydim(j) - ydim(i))
      end subroutine

 
