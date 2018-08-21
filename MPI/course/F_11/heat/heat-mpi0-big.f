        program heat
        include 'heat-mpi0-big.h' 

c additional variables for parallelization
        integer dims(1:2), coords(1:2)
        logical period(1:2)

c       naming: originally: i=0..imax,  now: i=is..ie,  icoord=0..idim-1
c       naming: originally: k=0..kmax,  now: k=ks..ke,  kcoord=0..kdim-1

        call MPI_INIT(ierror)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
        call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)

c test if size is 1, 4, 9, 16, 25, ... or 2, 8, 18, 32, 50, ...
        idim = int(sqrt(real(size)))
        kdim = int(sqrt(real(size)))
        if((idim*kdim).ne.(size)) then
          idim = int(sqrt(real(size/2)))
          kdim = 2*idim
        end if
        if((idim*kdim).ne.(size)) then
          if(rank.eq.0) then
            print *, 'Number of nodes must be n*n or 2*n*n with n>0'
          end if
          goto 500
        end if
        
        dims(1) = idim
        dims(2) = kdim
        period(1)=.false.
        period(2)=.false.

        call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,
     f                      period,.true.,comm,ierror)
        call MPI_COMM_RANK(comm, rank, ierror)
        call MPI_CART_COORDS(comm, rank, 2, coords, ierror) 
        icoord = coords(1) 
        kcoord = coords(2) 

c computing is:ie, ks:ke  ---  details see heat-mpi1-big.f 
        iinner1 = ((1 + imax - istart) - 2*b1 - 1) / idim + 1 
        in1 = (1 + imax - istart) - 2*b1 - idim * (iinner1 - 1) 
        if (icoord .lt. in1) then
          iinner = iinner1
          is = istart + icoord * iinner 
        else
          iinner = iinner1 - 1
          is = istart + in1 * iinner1 + (icoord-in1) * iinner 
        endif 
        iouter = iinner + 2*b1 
        ie = is + iouter - 1
c same for x coordinate: 
        kinner1 = ((1 + kmax - kstart) - 2*b1 - 1) / kdim + 1 
        kn1 = (1 + kmax - kstart) - 2*b1 - kdim * (kinner1 - 1) 
        if (kcoord .lt. kn1) then
          kinner = kinner1
          ks = kstart + kcoord * kinner 
        else
          kinner = kinner1 - 1
          ks = kstart + kn1 * kinner1 + (kcoord-kn1) * kinner 
        endif 
        kouter = kinner + 2*b1 
        ke = ks + kouter - 1

c The algorithm must be separated, because the array phi must no
c longer be a static array. It must be allocated dynamically. 

 
c It must be guaranteed that the array 'phin' is not an empty array: 
        if(((1 + imax - istart - 2*b1) .lt. idim) .or.
     f     ((1 + kmax - kstart - 2*b1) .lt. kdim)) then
          if(rank.eq.0) then
            print *, 'do not use more than',
     f               (1+imax-istart-2*b1)*(1+kmax-kstart-2*b1),
     f               'nodes to compute this application' 
          end if
          goto 500
        end if

        if (rank.eq.0) write (*,'(/a/a/)')
     f'!size iter-  wall clock time communication part  abort criterion'
     f, 
     f'!     ations    [seconds]    method [seconds]    meth. stride '//
     f'[seconds]'
 
         call algorithm(10,.false.)

c END of Programm
500     call MPI_FINALIZE(ierror)

c       ... no 'stop' statement, otherwise each process may issue
c           a 'stop' message on the output 
c       stop
        end
c
c------ SUBROUTINE ALGORITHM --------------------
c
        subroutine algorithm(stride,prt)
        include 'heat-mpi0-big.h' 
 
        integer stride
        logical prt
 
        double precision phi(is:ie,ks:ke), phin(is+b1:ie-b1,ks+b1:ke-b1)
        double precision dx,dy,dx2,dy2,dx2i,dy2i,dt,dphi,dphimax
        double precision dphimaxpartial
        double precision start_time, end_time
        double precision comm_time, criterion_time

        integer req(4), statuses(MPI_STATUS_SIZE,4)
        integer horizontal_border, vertical_border 

        call MPI_CART_SHIFT(comm, 0, 1, left, right, ierror)
        call MPI_CART_SHIFT(comm, 1, 1, lower,upper, ierror)

c create a MPI Vector
        call MPI_TYPE_VECTOR(kinner,b1,iouter, MPI_DOUBLE_PRECISION,
     f                       vertical_border,ierror)
        call MPI_TYPE_COMMIT(vertical_border, ierror)

        call MPI_TYPE_VECTOR(b1,iinner,iouter, MPI_DOUBLE_PRECISION,
     f                       horizontal_border,ierror)
        call MPI_TYPE_COMMIT(horizontal_border, ierror)
c         
        dx=1.d0/kmax
        dy=1.d0/imax
        dx2=dx*dx
        dy2=dy*dy
        dx2i=1.d0/dx2
        dy2i=1.d0/dy2
        dt=min(dx2,dy2)/4.d0
c start values 0.d0 
        do k=ks, min(ke,kmax-1)
        do i=max(1,is), min(ie,imax-1)
        phi(i,k)=0.d0
        enddo
        enddo
c start values 1.d0
       if (ke.eq.kmax) then 
        do i=is,ie
        phi(i,kmax)=1.d0
        enddo
       endif 
c start values dx
       if (is.eq.0) then 
c       phi(0,0)=0.d0
c       do k=1,kmax-1 
c         phi(0,k)=phi(0,k-1)+dx 
c         phi(imax,k)=phi(imax,k-1)+dx 
c       enddo 
c       ... substitute algorithmus by a code, 
c           that can be computed locally: 
        phi(0,ks) = ks*dx 
        do k=ks+1,min(ke,kmax-1)
         phi(0,k)=phi(0,k-1)+dx
        enddo
       endif 
       if (ie.eq.imax) then 
        phi(imax,ks) = ks*dx 
        do k=ks+1,min(ke,kmax-1)
         phi(imax,k)=phi(imax,k-1)+dx
        enddo
       endif 
 
c print array
       if (rank.eq.0) then
c       write (*,'(/,a)')
c    f  'Heat Conduction 2d'
c       write (*,'(/,4(a,1pg12.4))')
c    f  'dx =',dx,', dy =',dy,', dt =',dt,', eps =',eps
       endif 
        start_time = MPI_WTIME() 
        comm_time = 0
        criterion_time = 0 
 
c iteration
        do it=1,itmax
            dphimax=0.
            do k=ks+b1,ke-b1
             do i=is+b1,ie-b1
              dphi=(phi(i+1,k)+phi(i-1,k)-2.*phi(i,k))*dy2i
     f             +(phi(i,k+1)+phi(i,k-1)-2.*phi(i,k))*dx2i
              dphi=dphi*dt
              dphimax=max(dphimax,dphi)
              phin(i,k)=phi(i,k)+dphi
             enddo
            enddo
c save values
            do k=ks+b1,ke-b1
             do i=is+b1,ie-b1
              phi(i,k)=phin(i,k)
             enddo
            enddo
 
c for optimization: allreduce only each stride's loop:
          criterion_time = criterion_time - MPI_WTIME() 
            if (mod(it,stride) .eq. 0) then 
              if (size .gt. 1) then 
                dphimaxpartial = dphimax
                call MPI_ALLREDUCE(dphimaxpartial, dphimax, 1,
     f                        MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierror)
              endif 
              if(dphimax.lt.eps) goto 10
            endif 
          criterion_time = criterion_time + MPI_WTIME() 
 
          comm_time = comm_time - MPI_WTIME() 
c send and receive to/from upper/lower
            call MPI_IRECV(phi(is+b1,ks),1,horizontal_border,
     f                     lower,MPI_ANY_TAG,comm, req(1),ierror)
            call MPI_IRECV(phi(is+b1,ke),1,horizontal_border,
     f                     upper,MPI_ANY_TAG,comm, req(2),ierror)
            call MPI_ISEND(phi(is+b1,ke-b1),1,horizontal_border,
     f                     upper,0,comm, req(3),ierror)
            call MPI_ISEND(phi(is+b1,ks+b1),1,horizontal_border,
     f                     lower,0,comm, req(4),ierror)
            call MPI_WAITALL(4, req, statuses, ierror)

c send and receive to/from left/right
            call MPI_IRECV(phi(is,ks+b1),1,vertical_border,
     f                     left,MPI_ANY_TAG,comm, req(1),ierror)
            call MPI_IRECV(phi(ie,ks+b1),1,vertical_border,
     f                     right,MPI_ANY_TAG,comm, req(2),ierror)
            call MPI_ISEND(phi(ie-b1,ks+b1),1,vertical_border,
     f                     right,0,comm, req(3),ierror)
            call MPI_ISEND(phi(is+b1,ks+b1),1,vertical_border,
     f                     left,0,comm, req(4),ierror)
            call MPI_WAITALL(4, req, statuses, ierror)
          comm_time = comm_time + MPI_WTIME() 
 
        enddo
10      continue
        criterion_time = criterion_time + MPI_WTIME() 
        end_time = MPI_WTIME() 
c print array
        if (prt) call heatpr(phi) 
       if (rank.eq.0) then 
         write(*,
     f   '(''!''i4,1x,i6,2x,g12.4,4x,i2,2x,g12.4,2x,i2,2x,i6,1x,g12.4)')
     f    size, it, end_time - start_time, 
     f    1, comm_time,
     f    1, stride, criterion_time 
       endif 
c
        return
        end
c
c --------------------------------------------------------------- 
c
        subroutine heatpr(phi)
        include 'heat-mpi0-big.h' 
        double precision phi(is:ie,ks:ke)
        do ic = 0, idim-1 
          i_first = is + b1
          if (ic.eq.0) i_first = i_first - b1
          i_last = ie - b1
          if (ic.eq.idim-1) i_last = i_last + b1
          do kc = 0, kdim-1
            if ((ic.eq.icoord) .and. (kc.eq.kcoord)) then 
              if (kc.eq.0) write(*,'(/a,i3,a)')
     f                  'printing the ',ic,'th horizontal block'
              k_first = ks + b1
              if (kc.eq.0) k_first = k_first - b1
              k_last = ke - b1
              if (kc.eq.kdim-1) k_last = k_last + b1
              do k=k_first,k_last
                write(*,'(i3,'':'',1000(1x,f4.2))')
     f                    k, (phi(i,k),i=i_first,i_last) 
              enddo
            endif 
            call MPI_BARRIER(comm,ierror)
          enddo 
        enddo 
        return
        end
