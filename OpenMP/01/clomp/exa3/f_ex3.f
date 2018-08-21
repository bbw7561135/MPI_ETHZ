      module nmod
      integer niter
      parameter (m=3,n=6,p=2)
      end module 

      use nmod
      common /matrices/ c(m,p),a(m,n),b(n,p) 
      integer a,b,c,i,j

      niter = 3

!initialize the source matrices
      do i=1, n
	  do j=1, m
	     a(j,i) = i+j 
	  enddo
      enddo

      do i=1, p
	  do j=1, n
	     b(j,i) = j - (i*2)
	  enddo
      enddo
      
      call run(niter - 1)
      end

      subroutine run(argexpr)
      use nmod
      integer argexpr
      common /matrices/ c(m,p),a(m,n),b(n,p) 
      integer a,b,c,i,j,k,l

!$omp parallel private(i,j,k,l) shared(c,a,b,argexpr,niter)
      do i=1, niter

!reset
!$omp do
         do j=1, p
	     do k=1, m
	         c(k, j) = 0
	     enddo
	 enddo

!multiply
!$omp do
         do j=1, p
	     do k=1, m
	         do l=1, n
		     c(k, j) = c(k, j) + (argexpr * a(k, l) * b(l, j))
	         enddo
	     enddo
	 enddo

!$omp single
      print *, "iter ", i
      print *, ""
      print *, "array ", c 
      print *, ""
!$omp end single
      enddo
!$omp end parallel



      end subroutine


