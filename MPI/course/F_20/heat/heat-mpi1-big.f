        program heat
        include 'heat-mpi1-big.h' 

c additional variables for parallelization
        integer dims(1:2), coords(1:2)
        logical period(1:2)
        integer communication_method, end_method, stride
        logical prt

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

c whole y indecees   |------------- isize = 1+imax-istart ------------|
c    start/end index  ^-istart                                  imax-^
c 
c 1. interval        |--- iouter1---|   
c                    |--|--------|--|
c                     b1  iinner1 b1 
c    start/end index  ^-is      ie-^ 
c 2. interval                 |--|--------|--|   
c 3. interval                          |--|--------|--| 
c 4. interval                                   |--|-------|--| 
c                                                   iinner0 
c 5. interval = idim's interval                         |--|-------|--|
c
c In each iteration on each interval, the inner area is computed
c by using the values of the last iteration in the whole outer area. 
c
c icoord = number of the interval - 1
c 
c To fit exactly into isize, we use in1 intervals of with iinner1
c and (idim-in1) intervals of with (iinner1 - 1) 
c
c isize <= 2*b1 + idim * iinner1  
c
c ==>   iinner1 >= (isize - 2*b1) / idim
c ==>   smallest valid integer value:

        iinner1 = ((1 + imax - istart) - 2*b1 - 1) / idim + 1 

c isize - 2*b1 = in1 * iinner1 + (idim - in1) * (iinner1 - 1)
c ==>
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

c       print *, 'rank=', rank, 'icoord=',icoord,'is,ie=',is,ie,
c    f                          'iouter=',iouter,
c    f                          'kcoord=',kcoord,'ks,ke=',ks,ke, 
c    f                          'kouter=',kouter
 
        if (rank.eq.0) write (*,'(/a/a/)')
     f'!size iter-  wall clock time communication part  abort criterion'
     f, 
     f'!     ations    [seconds]    method [seconds]    meth. stride '//
     f'[seconds]'
 
c       if (.true.) then 
        if (.false.) then 
         if (rank.eq.0) then
          write (*,*)'type communication_method [1..7] end_method [1-3]'
     f               //' stride [>0] print_flag [t,f]' 
          read  (*,*) communication_method, end_method, stride, prt  
         endif 
         CALL MPI_BCAST(communication_method,1,MPI_INTEGER,
     f                  0,comm,ierror) 
         CALL MPI_BCAST(end_method,1,MPI_INTEGER,0,comm,ierror) 
         CALL MPI_BCAST(stride,1,MPI_INTEGER,0,comm,ierror) 
         CALL MPI_BCAST(prt,1,MPI_LOGICAL,0,comm,ierror) 
         call algorithm(communication_method, end_method, stride, prt) 
        else 
         prt = .false. 
c        prt = .true.
         if (size.gt.1) then 
          call algorithm(1,1,179,prt)
          call algorithm(2,1,179,prt)
          call algorithm(3,1,179,prt)
          call algorithm(4,1,179,prt)
          call algorithm(5,1,179,prt)
          call algorithm(6,1,179,prt)
          call algorithm(7,1,179,prt)
 
          call algorithm(5,1,1,prt)
          call algorithm(5,1,10,prt)
          call algorithm(5,1,20,prt)
          call algorithm(5,1,80,prt)
          call algorithm(5,1,98,prt)
          call algorithm(5,1,179,prt)
          call algorithm(5,2,0,prt)
          call algorithm(5,3,0,prt)
 
          call algorithm(7,1,  1,prt)
          call algorithm(7,1, 80,prt)
          call algorithm(7,1, 98,prt)
          call algorithm(7,1,179,prt)
          call algorithm(7,1, 80,prt)
          call algorithm(7,1, 98,prt)
          call algorithm(7,1,179,prt)
         else
          call algorithm(0,1,1,prt) 
          call algorithm(0,1,80,prt) 
          call algorithm(0,1,98,prt) 
          call algorithm(0,1,179,prt) 
         endif 
        endif 

c END of Programm
500     call MPI_FINALIZE(ierror)

c       ... no 'stop' statement, otherwise each process may issue
c           a 'stop' message on the output 
c       stop
        end
c
c
c------ SUBROUTINE ALGORITHM --------------------
c
        subroutine algorithm(communication_method,end_method,stride,prt)
        include 'heat-mpi1-big.h' 
 
        integer communication_method, end_method, stride
        logical prt
 
        double precision phi(is:ie,ks:ke), phin(is+b1:ie-b1,ks+b1:ke-b1)
        double precision dx,dy,dx2,dy2,dx2i,dy2i,dt,dphi,dphimax
        double precision dphimaxpartial
        double precision start_time, end_time
        double precision comm_time, criterion_time

        integer done_local, reached_local
        integer it_stopping 

c only for communication7: 
        integer sendcounts(0:size-1), sdispls(0:size-1) 
        integer recvcounts(0:size-1), rdispls(0:size-1) 
        double precision phiright(ks+1:ke-1)
        double precision phileft(ks+1:ke-1)
        double precision phirecvright(ks+1:ke-1)
        double precision phirecvleft(ks+1:ke-1)
        integer dp_size 
c ------- 
 
        call MPI_CART_SHIFT(comm, 0, 1, left, right, ierror)
        call MPI_CART_SHIFT(comm, 1, 1, lower,upper, ierror)

c create a MPI Vector
        call MPI_TYPE_VECTOR(kinner,b1,iouter, MPI_DOUBLE_PRECISION,
     f                       vertical_border,ierror)
        call MPI_TYPE_COMMIT(vertical_border, ierror)

        call MPI_TYPE_VECTOR(b1,iinner,iouter, MPI_DOUBLE_PRECISION,
     f                       horizontal_border,ierror)
        call MPI_TYPE_COMMIT(horizontal_border, ierror)
 
 
c for communication7 with MPI_ALLTOALLV 
        call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, dp_size, ierror)
        do i=0, size-1
          sendcounts(i) = 0
          sdispls(i) = 0 
          recvcounts(i) = 0
          rdispls(i) = 0 
        enddo 
 
        if (upper.ne.MPI_PROC_NULL) then 
          sendcounts(upper) = iinner
          call MPI_ADDRESS(phi(is+1,ke-1), sdispls(upper), ierror)
          sdispls(upper) = sdispls(upper) / dp_size 
          recvcounts(upper) = iinner
          call MPI_ADDRESS(phi(is+1,ke)  , rdispls(upper), ierror)
          rdispls(upper) = rdispls(upper) / dp_size 
        endif 
        if (lower.ne.MPI_PROC_NULL) then 
          sendcounts(lower) = iinner
          call MPI_ADDRESS(phi(is+1,ks+1), sdispls(lower), ierror)
          sdispls(lower) = sdispls(lower) / dp_size 
          recvcounts(lower) = iinner
          call MPI_ADDRESS(phi(is+1,ks)  , rdispls(lower), ierror)
          rdispls(lower) = rdispls(lower) / dp_size 
        endif 
        if (right.ne.MPI_PROC_NULL) then 
          sendcounts(right) = kinner
          call MPI_ADDRESS(phiright(ks+1), sdispls(right), ierror)
          sdispls(right) = sdispls(right) / dp_size 
          recvcounts(right) = kinner
          call MPI_ADDRESS(phirecvright(ks+1), rdispls(right), ierror)
          rdispls(right) = rdispls(right) / dp_size 
        endif 
        if (left.ne.MPI_PROC_NULL) then 
          sendcounts(left) = kinner
          call MPI_ADDRESS(phileft(ks+1), sdispls(left), ierror)
          sdispls(left) = sdispls(left) / dp_size 
          recvcounts(left) = kinner
          call MPI_ADDRESS(phirecvleft(ks+1), rdispls(left), ierror)
          rdispls(left) = rdispls(left) / dp_size 
        endif 
c ----------- 
         
        done_local = 0
        lower_left_tag = 0
        upper_right_tag = 0 
        if (end_method.eq.2) then 
          reached_local = 0
        else if(end_method.eq.3) then 
          it_stopping = 0 
          sts_upper(MPI_TAG) = 0 
          sts_lower(MPI_TAG) = 0 
          sts_right(MPI_TAG) = 0 
          sts_left (MPI_TAG) = 0 
        endif 
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
          if (done_local.eq.0) then 
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
          endif 
c for optimization: allreduce only each stride's loop:
          criterion_time = criterion_time - MPI_WTIME() 
          goto (101, 102, 103), end_method 
101        continue 
c          ... end_method == 1 
            if (mod(it,stride) .eq. 0) then 
              if (size .gt. 1) then 
                dphimaxpartial = dphimax
                call MPI_ALLREDUCE(dphimaxpartial, dphimax, 1,
     f                        MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierror)
              endif 
              if(dphimax.lt.eps) goto 110
            endif 
           goto 199 
 
102        continue 
c          ... end_method == 2 
            if (done_local.eq.1) then
c             ... After sending "done_local" in the last iteration: 
c                 if a neighbor is also in the state 'done_local' then 
c                 the neighbor now knows that this node is also in
c                 'done_local' and both can stop to communicate:
              if (sts_upper(MPI_TAG).eq.1) upper = MPI_PROC_NULL
              if (sts_lower(MPI_TAG).eq.1) lower = MPI_PROC_NULL
              if (sts_right(MPI_TAG).eq.1) right = MPI_PROC_NULL
              if (sts_left (MPI_TAG).eq.1) left  = MPI_PROC_NULL
c             ... if there is no more any local computation
c                 nor communication then the process can be stopped: 
              if ((upper.eq.MPI_PROC_NULL).and.(lower.eq.MPI_PROC_NULL)
     f        .and.(right.eq.MPI_PROC_NULL).and.(left.eq.MPI_PROC_NULL))
     f          goto 110 
            endif 
c           ... "reached_local" means that the heat wave has reached
c               the local area of phi.
c               To test "reached_local" is necessary to prohibit that
c               a node stops before it has computed anything due
c               to its starting values that are all zero. 
            if (dphimax.gt.eps) reached_local = 1 
            if ((reached_local.eq.1).and.(dphimax.lt.eps)) then
              done_local = 1
              lower_left_tag = 1
              upper_right_tag = 1 
            endif 
           goto 199 
 
103        continue 
c          ... end_method == 3 
            upper_right_tag = 0 
c           ... outgoing from node (0,0) the local_done criterion
c               is propageded with AND operation to the upper-right node
            if ((dphimax.lt.eps) .and.
     f       ((lower.eq.MPI_PROC_NULL).or.(sts_lower(MPI_TAG).eq.1))
     f       .and.((left.eq.MPI_PROC_NULL).or.(sts_left(MPI_TAG).eq.1)))
     f         upper_right_tag = 1 
c           ... after the local_done criterion has reached at the
c               uppermost-rightmost node in the interval of the 
c               last (idim+kdim) iterations, it is sent back to
c               the node (0,0) to inform all nodes about the global stop
            if ((icoord.eq.idim-1) .and. (kcoord.eq.kdim-1)) then
             if ((dphimax.lt.eps) .and.
     f       ((lower.eq.MPI_PROC_NULL).or.(sts_lower(MPI_TAG).eq.1))
     f       .and.((left.eq.MPI_PROC_NULL).or.(sts_left(MPI_TAG).eq.1)))
     f         lower_left_tag = 1 
            else 
             if (((upper.eq.MPI_PROC_NULL).or.(sts_upper(MPI_TAG).eq.1))
     f       .and.
     f       ((right.eq.MPI_PROC_NULL).or.(sts_right(MPI_TAG).eq.1)))
     f         lower_left_tag = 1 
            endif 
            it_stopping = it_stopping + lower_left_tag
c           ... then the global stop is done synchronously after all
c               nodes are informed: 
            if (it_stopping .eq. (icoord+kcoord+1)) goto 110 
           goto 199

110       continue 
           criterion_time = criterion_time + MPI_WTIME() 
          goto 10 
199       continue 
           criterion_time = criterion_time + MPI_WTIME() 
 
c         ... The following "if" statement is not necessary. 
c             It optimizes the case "size=1":
          if (size .gt. 1) then 
           comm_time = comm_time - MPI_WTIME() 
           goto (201,202,203,204,205,206,207), communication_method 
201         call communication1(phi) 
            goto 299
202         call communication2(phi) 
            goto 299
203         call communication3(phi) 
            goto 299
204         call communication4(phi) 
            goto 299
205         call communication5(phi) 
            goto 299
206         call communication6(phi) 
            goto 299
207         continue 
             phileft(ks+1:ke-1)=phi(is+1,ks+1:ke-1)
             phiright(ks+1:ke-1)=phi(ie-1,ks+1:ke-1)
c            ... to prevent compiler optimizations because the compiler
c                does not see the input arguments of ALLTOALLV:
             call MPI_ADDRESS(phi,iaddr_dummy,ierror) 
             call MPI_ADDRESS(phileft,iaddr_dummy,ierror) 
             call MPI_ADDRESS(phiright,iaddr_dummy,ierror) 
             call MPI_ALLTOALLV(
     f           MPI_BOTTOM, sendcounts, sdispls, MPI_DOUBLE_PRECISION, 
     f           MPI_BOTTOM, recvcounts, rdispls, MPI_DOUBLE_PRECISION, 
     f           comm, ierror) 
c            ... to prevent compiler optimizations because the compiler
c                does not see the output arguments of ALLTOALLV:
             call MPI_ADDRESS(phi,iaddr_dummy,ierror) 
             call MPI_ADDRESS(phirecvleft,iaddr_dummy,ierror) 
             call MPI_ADDRESS(phirecvright,iaddr_dummy,ierror) 
             if (left.ne.MPI_PROC_NULL)
     f         phi(is,ks+1:ke-1)=phirecvleft(ks+1:ke-1) 
             if (right.ne.MPI_PROC_NULL) 
     f         phi(ie,ks+1:ke-1)=phirecvright(ks+1:ke-1) 
            goto 299
299        continue 
          comm_time = comm_time + MPI_WTIME() 
          endif 
        enddo
10      continue
        end_time = MPI_WTIME() 
c print array
        if (prt) call heatpr(phi) 
       if (rank.eq.0) then 
         write(*,
     f   '(''!''i4,1x,i6,2x,g12.4,4x,i2,2x,g12.4,2x,i2,2x,i6,1x,g12.4)')
     f    size, it, end_time - start_time, 
     f    communication_method, comm_time,
     f    end_method, stride, criterion_time 
       endif 
c
        return
        end
c
c
c
        subroutine heatpr(phi)
        include 'heat-mpi1-big.h' 
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
c
c
c
c------ different communication routines
c
c------ SUBROUTINE COMMUNICATION1 -----------------------
c
        subroutine communication1(phi)
        include 'heat-mpi1-big.h' 
        double precision phi(is:ie,ks:ke)

        integer req(4)

c send and receive to/from upper/lower
c       ... The following "if" statement is not necessary. 
c           It optimizes the case "kdim=1":
        if (kdim.gt.1) then 
         call MPI_IRECV(phi(is+b1,ks),1,horizontal_border,
     f                  lower,MPI_ANY_TAG,comm, req(2),ierror)
         call MPI_IRECV(phi(is+b1,ke),1,horizontal_border,
     f                  upper,MPI_ANY_TAG,comm, req(1),ierror)
         call MPI_ISEND(phi(is+b1,ke-b1),1,horizontal_border,
     f                  upper,upper_right_tag,comm, req(3),ierror)
         call MPI_ISEND(phi(is+b1,ks+b1),1,horizontal_border,
     f                  lower,lower_left_tag,comm, req(4),ierror)
         call MPI_WAITALL(4, req, sts_upper, ierror)
        endif 

c send and receive to/from left/right
c       ... The following "if" statement is not necessary. 
c           It optimizes the case "idim=1":
        if (idim.gt.1) then 
         call MPI_IRECV(phi(is,ks+b1),1,vertical_border,
     f                  left,MPI_ANY_TAG,comm, req(2),ierror)
         call MPI_IRECV(phi(ie,ks+b1),1,vertical_border,
     f                  right,MPI_ANY_TAG,comm, req(1),ierror)
         call MPI_ISEND(phi(ie-b1,ks+b1),1,vertical_border,
     f                  right,upper_right_tag,comm, req(3),ierror)
         call MPI_ISEND(phi(is+b1,ks+b1),1,vertical_border,
     f                  left,lower_left_tag,comm, req(4),ierror)
         call MPI_WAITALL(4, req, sts_right, ierror)
        endif 

        return
        end

c
c
c------ SUBROUTINE COMMUNICATION2 -----------------------
c
        subroutine communication2(phi)
        include 'heat-mpi1-big.h' 
        double precision phi(is:ie,ks:ke)

c CAUTION: This routine requires that b1 = 1 !!! 

        integer req(4)

c send and receive to/from upper/lower
c       ... The following "if" statement is not necessary. 
c           It optimizes the case "kdim=1":
        if (kdim.gt.1) then 
         call MPI_IRECV(phi(is+1,ks),iinner,MPI_DOUBLE_PRECISION,
     f                  lower,MPI_ANY_TAG,comm, req(2),ierror)
         call MPI_IRECV(phi(is+1,ke),iinner,MPI_DOUBLE_PRECISION,
     f                  upper,MPI_ANY_TAG,comm, req(1),ierror)
         call MPI_ISEND(phi(is+1,ke-1),iinner,MPI_DOUBLE_PRECISION,
     f                  upper,upper_right_tag,comm, req(3),ierror)
         call MPI_ISEND(phi(is+1,ks+1),iinner,MPI_DOUBLE_PRECISION,
     f                  lower,lower_left_tag,comm, req(4),ierror)
         call MPI_WAITALL(4, req, sts_upper, ierror)
        endif 

c send and receive to/from left/right
c       ... The following "if" statement is not necessary. 
c           It optimizes the case "idim=1":
        if (idim.gt.1) then 
         call MPI_IRECV(phi(is,ks+1),1,vertical_border,
     f                  left,MPI_ANY_TAG,comm, req(2),ierror)
         call MPI_IRECV(phi(ie,ks+1),1,vertical_border,
     f                  right,MPI_ANY_TAG,comm, req(1),ierror)
         call MPI_ISEND(phi(ie-1,ks+1),1,vertical_border,
     f                  right,upper_right_tag,comm, req(3),ierror)
         call MPI_ISEND(phi(is+1,ks+1),1,vertical_border,
     f                  left,lower_left_tag,comm, req(4),ierror)
         call MPI_WAITALL(4, req, sts_right, ierror)
        endif 

        return
        end

c
c
c------ SUBROUTINE COMMUNICATION3 -----------------------
c
        subroutine communication3(phi)
        include 'heat-mpi1-big.h' 
        double precision phi(is:ie,ks:ke)

c CAUTION: This routine requires that b1 = 1 !!! 

        integer req(8)
        double precision phiright(ks+1:ke-1)
        double precision phileft(ks+1:ke-1)

c       ... If "WAITALL(8,...)" is used then the following assignments
c           must be done before the upper and lower "ISEND", because
c           the corners phi(is+1,ks+1), phi(is+1,ke-1), phi(ie-1,ks+1),
c           and phi(ie-1,ke-1) are sent upper/lower and left/right
c           and because the MPI 1.1 has forbidden to "ISEND" twice the
c           same memory location, see MPI 1.1, page 40, lines 44+45.
        phileft(ks+1:ke-1)=phi(is+1,ks+1:ke-1)
        phiright(ks+1:ke-1)=phi(ie-1,ks+1:ke-1)

c send and receive to/from upper/lower
        call MPI_IRECV(phi(is+1,ks),iinner,MPI_DOUBLE_PRECISION,
     f                 lower,MPI_ANY_TAG,comm, req(2),ierror)
        call MPI_IRECV(phi(is+1,ke),iinner,MPI_DOUBLE_PRECISION,
     f                 upper,MPI_ANY_TAG,comm, req(1),ierror)
        call MPI_ISEND(phi(is+1,ke-1),iinner,MPI_DOUBLE_PRECISION,
     f                 upper,upper_right_tag,comm, req(5),ierror)
        call MPI_ISEND(phi(is+1,ks+1),iinner,MPI_DOUBLE_PRECISION,
     f                 lower,lower_left_tag,comm, req(6),ierror)

c send and receive to/from left/right
        call MPI_IRECV(phi(is,ks+1),1,vertical_border,
     f                 left,MPI_ANY_TAG,comm, req(4),ierror)
        call MPI_IRECV(phi(ie,ks+1),1,vertical_border,
     f                 right,MPI_ANY_TAG,comm, req(3),ierror)
        call MPI_ISEND(phiright(ks+1),kinner,MPI_DOUBLE_PRECISION,
     f                 right,upper_right_tag,comm, req(7),ierror)
        call MPI_ISEND(phileft(ks+1),kinner,MPI_DOUBLE_PRECISION,
     f                 left,lower_left_tag,comm, req(8),ierror)

c waiting till all send/receive finished
        call MPI_WAITALL(8, req, sts_array, ierror)

        return
        end

c
c
c------ SUBROUTINE COMMUNICATION4 -----------------------
c
        subroutine communication4(phi)
        include 'heat-mpi1-big.h' 
        double precision phi(is:ie,ks:ke)

c CAUTION: This routine requires that b1 = 1 !!! 

        integer req(8)
        double precision phiright(ks+1:ke-1)
        double precision phileft(ks+1:ke-1)

c       ... If "WAITALL(8,...)" is used then the following assignments
c           must be done before the upper and lower "ISEND", because
c           the corners phi(is+1,ks+1), phi(is+1,ke-1), phi(ie-1,ks+1),
c           and phi(ie-1,ke-1) are sent upper/lower and left/right
c           and because the MPI 1.1 has forbidden to "ISEND" twice the
c           same memory location, see MPI 1.1, page 40, lines 44+45.
        phileft(ks+1:ke-1)=phi(is+1,ks+1:ke-1)
        phiright(ks+1:ke-1)=phi(ie-1,ks+1:ke-1)

c send and receive to/from upper/lower
        call MPI_IRECV(phi(is+1,ks),iinner,MPI_DOUBLE_PRECISION,
     f                 lower,MPI_ANY_TAG,comm, req(2),ierror)
        call MPI_ISEND(phi(is+1,ke-1),iinner,MPI_DOUBLE_PRECISION,
     f                 upper,upper_right_tag,comm, req(5),ierror)
        call MPI_IRECV(phi(is+1,ke),iinner,MPI_DOUBLE_PRECISION,
     f                 upper,MPI_ANY_TAG,comm, req(1),ierror)
        call MPI_ISEND(phi(is+1,ks+1),iinner,MPI_DOUBLE_PRECISION,
     f                 lower,lower_left_tag,comm, req(6),ierror)

c send and receive to/from left/right
        call MPI_IRECV(phi(is,ks+1),1,vertical_border,
     f                 left,MPI_ANY_TAG,comm, req(4),ierror)
        call MPI_ISEND(phiright(ks+1),kinner,MPI_DOUBLE_PRECISION,
     f                 right,upper_right_tag,comm, req(7),ierror)
        call MPI_IRECV(phi(ie,ks+1),1,vertical_border,
     f                 right,MPI_ANY_TAG,comm, req(3),ierror)
        call MPI_ISEND(phileft(ks+1),kinner,MPI_DOUBLE_PRECISION,
     f                 left,lower_left_tag,comm, req(8),ierror)

c waiting till all send/receive finished
        call MPI_WAITALL(8, req, sts_array, ierror)

        return
        end
c
c
c------ SUBROUTINE COMMUNICATION5 -----------------------
c
        subroutine communication5(phi)
        include 'heat-mpi1-big.h' 
        double precision phi(is:ie,ks:ke)

c CAUTION: This routine requires that b1 = 1 !!! 

        double precision phiright(ks+1:ke-1)
        double precision phileft(ks+1:ke-1)

c send and receive to/from upper/lower
c       ... The following "if" statement is not necessary. 
c           It optimizes the case "kdim=1":
        if (kdim.gt.1) then 
         call MPI_SENDRECV(phi(is+1,ke-1),iinner,
     f                     MPI_DOUBLE_PRECISION,upper,upper_right_tag,
     f                     phi(is+1,ks),iinner,
     f                     MPI_DOUBLE_PRECISION,lower,MPI_ANY_TAG,
     f                     comm,sts_lower,ierror)
         call MPI_SENDRECV(phi(is+1,ks+1),iinner,
     f                     MPI_DOUBLE_PRECISION,lower,lower_left_tag,
     f                     phi(is+1,ke),iinner,
     f                     MPI_DOUBLE_PRECISION,upper,MPI_ANY_TAG,
     f                     comm,sts_upper,ierror)
        endif 

c send and receive to/from left/right
c       ... The following "if" statement is not necessary. 
c           It optimizes the case "idim=1":
        if (idim.gt.1) then 
         phileft(ks+1:ke-1)=phi(is+1,ks+1:ke-1)
         phiright(ks+1:ke-1)=phi(ie-1,ks+1:ke-1)
         call MPI_SENDRECV(phiright(ks+1),kinner,
     f                     MPI_DOUBLE_PRECISION,right,upper_right_tag,
     f                     phi(is,ks+1),1,
     f                     vertical_border,left,MPI_ANY_TAG,
     f                     comm,sts_left,ierror)
         call MPI_SENDRECV(phileft(ks+1),kinner,
     f                     MPI_DOUBLE_PRECISION,left,lower_left_tag,
     f                     phi(ie,ks+1),1,
     f                     vertical_border,right,MPI_ANY_TAG,
     f                     comm,sts_right,ierror)
        endif 

        return
        end
c
c
c------ SUBROUTINE COMMUNICATION6 -----------------------
c
        subroutine communication6(phi)
        include 'heat-mpi1-big.h' 
        double precision phi(is:ie,ks:ke)

c CAUTION: This routine requires that b1 = 1 !!! 

        double precision phiright(ks+1:ke-1)
        double precision phileft(ks+1:ke-1)
        double precision phirecvright(ks+1:ke-1)
        double precision phirecvleft(ks+1:ke-1)

c send and receive to/from upper/lower
c       ... The following "if" statement is not necessary. 
c           It optimizes the case "kdim=1":
        if (kdim.gt.1) then 
         call MPI_SENDRECV(phi(is+1,ke-1),iinner,
     f                     MPI_DOUBLE_PRECISION,upper,upper_right_tag,
     f                     phi(is+1,ks),iinner,
     f                     MPI_DOUBLE_PRECISION,lower,MPI_ANY_TAG,
     f                     comm,sts_lower,ierror)
         call MPI_SENDRECV(phi(is+1,ks+1),iinner,
     f                     MPI_DOUBLE_PRECISION,lower,lower_left_tag,
     f                     phi(is+1,ke),iinner,
     f                     MPI_DOUBLE_PRECISION,upper,MPI_ANY_TAG,
     f                     comm,sts_upper,ierror)
        endif 

c send and receive to/from left/right
c       ... The following "if" statement is not necessary. 
c           It optimizes the case "idim=1":
        if (idim.gt.1) then 
         phileft(ks+1:ke-1)=phi(is+1,ks+1:ke-1)
         phiright(ks+1:ke-1)=phi(ie-1,ks+1:ke-1)
         call MPI_SENDRECV(phiright(ks+1),kinner,
     f                     MPI_DOUBLE_PRECISION,right,upper_right_tag,
     f                     phirecvleft(ks+1),kinner,
     f                     MPI_DOUBLE_PRECISION,left,MPI_ANY_TAG,
     f                     comm,sts_left,ierror)
         call MPI_SENDRECV(phileft(ks+1),kinner,
     f                     MPI_DOUBLE_PRECISION,left,lower_left_tag,
     f                     phirecvright(ks+1),kinner,
     f                     MPI_DOUBLE_PRECISION,right,MPI_ANY_TAG,
     f                     comm,sts_right,ierror)
         if (left.ne.MPI_PROC_NULL)
     f     phi(is,ks+1:ke-1)=phirecvleft(ks+1:ke-1) 
         if (right.ne.MPI_PROC_NULL) 
     f     phi(ie,ks+1:ke-1)=phirecvright(ks+1:ke-1) 
        endif 

        return
        end

 
!====================================================================== 
 
! Comments:
! 
! * static arrays must be substituted by dynamically allocated arrays,
!   therefore additional subroutine "algorithm" that allocates "phi"
! 
! * block partitioning: communication time = 
!                         4 * sqrt(area) / bandwidth / sqrt(#PE)
!                         + 4 * latency              ===========
!   versus
!   horizontal stripes: communication time = 
!                         2 * sqrt(area) / bandwidth
!                         + 2 * latency
!                         ===
! 
! * MPI topology functions help to organise the indicees and coordinates
! 
! * communication choices:
! 
!   - goals
!     -- message passing
!     -- no deadlocks
!     -- correct MPI code
!     -- efficient
! 
!   - methods
!     -- non-blocking point-to-point  (examples 1-4)
!     -- MPI_SENDRECV  (5-6)
!     -- MPI_ALLTOALLV (7)
! 
!   - details and other choices
!     -- strided datatype (1-5) versus local copying (3-7) 
!     -- separate vertical and horizontal sending to
!        prohibit sinding of the corners into two directions
!        simultanously (see MPI 1.1, , page 40, lines 44+45) (1+2)
!        versus usage of local copies (3+4) 
!     -- IRECV in both directions (of one dimension) and afterwards
!        the two ISENDs (1,2,3)
!        versus two pairs of IRECV+ISEND (4)
! 
! * abort criterion
! 
!   - original (sequential) algorithm
! 
!        do 
!          dphimax=0.
!          do k=1,kmax-1
!            do i=1,imax-1
!              dphi=...
!              dphimax=max(dphimax,dphi)
!            enddo
!          enddo
!          ... 
!          if(dphimax.lt.eps) goto 10
!        enddo
!     10 continue 
! 
!   - goals
!     -- efficient
!     -- similar to original algorithm
!     -- numerically and MPI-technically correct 
! 
!   - solutions
! 
!     1a. With MPI_ALLREDUCE (identical to sequential algorithm),
!         i.e. stride == 1.
!                     
!         DISADVANTAGE: 
!         - Although MPI_ALLREDUCE is implemented mostly
!           with a binary tree algorithm, it is very expensive. 
!  
!     1b. With MPI_ALLREDUCE, but not in each iteration, i.e. stride>1.
!                     
!         DISADVANTAGE: 
!         - Between zero and (stride-1) additional iterations must be
!           computed. 
! 
!     2.  Stop computation if criterion is reached locally, and
!         stop communication with a neighbor if the neighbor
!         has also reached the criterion, and
!         stop all work on a process if computation locally 
!         and communication to all neighbors is stopped.
! 
!         The information "criterion is reached locally"
!         is sent to the neighbors via the MPI tag.
!  
!         DISADVANTAGES: 
!         - If the criterion is reached locally
!           in one process before e.g. a wave has reached this
!           node then this node will stop too early.
!           Workaround for this problem:
!           Testing the boart criterion is delayed until
!           the first time the processor has done significant work,
!           i.e. until dphimax is greater than eps the first time.
!         - The processes are not stopped simultanously at the same
!           iteration.   
!         - This method cannot be combined with the communication
!           method "MPI_ALLTOALLV" because MPI tags are used. 
! 
!     3.  A two-phase protocol:
! 
!         Phase one: The "locally done" criterion is propagated
!                    from the (0,0)-node to the rightmost-uppermost 
!                    node by using the MPI tag in the right and 
!                    upper direction.
!                    Each node sends "true" only if it received "true"
!                    and if the criterion is reached locally.
!                    If the criterion is losed again, then the
!                    node sends "false" again.
!                    This phase ends, if the rightmost-uppermost
!                    node has reached the criterion locally and 
!                    has received "true" from left and lower nodes.
!         Phase two: As soon as the rightmost-uppermost has finished
!                    phase one, this node sends into the directions
!                    left and lower the MPI tag "true".
!                    Each node propagates the information into the
!                    same direction.
!                    This phase is finished on each node simultanously
!                    (icoord+kcoord) steps after the tag has reached
!                    the node.
!         After phase two all processes can stop simultanously.
!                     
!         DISADVANTAGES: 
!         - This method cannot be combined with the communication
!           method "MPI_ALLTOALLV" because MPI tags are used. 
!         - In the worst case 2*(idim+kdim-2) additional
!           iterations are computed, in the best case 
!           1*(idim+kdim-2) additional iterations. 
! 
! * Best choices are machine dependent:
! 
!   - CRAY T3E:        ALLTOALLV      with  ALLREDUCE, stride=98 
!  
!   - Hitachi SR2201:  non-blocking   with  end_method 2 or 3
!                      or SENDRECV  
!
! * References to MPI chapters and parallelization technics:
! 
!   - blocking MPI_SENDRECV       --> communication5 + communication6 
!   - non-blocking point-to-point --> communication1 - communication4
!   - collective operations      
!      -- MPI_ALLTOALLV           --> communication_method 7 
!                                     lines 447-468 (label 207)
!                                     l. 201-208, 224-265 (preparation) 
!      -- MPI_BCAST               --> lines 133-137
!      -- MPI_ALLREDUCE           --> end_method 1 (line 357)
!   - derived datatypes           --> lines 215-221 
!   - usage of MPI_BOTTOM         --> see MPI_ALLTOALLV references
!     and MPI_ADDRESS 
!   - problems with Fortran       --> lines 452-455, 461-463, solved
!     bindings for MPI                with dummy calls to MPI_ADDRESS
!     (MPI-2, pages 289f)             (used as DD, MPI-2 p.290, l.8-32),
!                                     necessary with CRAY T3E f90-comp. 
!   - block data decomposition    --> lines 44-101, 192
!   - MPI topologies              --> lines 32-42, 211-212, and
!                                     "comm" is used in each
!                                     MPI communication routine
!   - MPI init+finalize           --> lines 14-16, 176
!   - MPI_PROC_NULL               --> assignment:
!     (MPI 1.1 Chapter. 3.11,         - lines 211-212 (MPI_CART_SHIFT)
!     pages 60f)                      - lines 371-374 (explicit)
!                                     usage:
!                                     - all point-to-point
!                                       communication routines
!                                     - explicit tests:
!                                       lines 234-265, 377+378, 400-415,
!                                             464+466, 802+804 
!   - timer MPI_WTIME             --> l. 326+474, 433+470, 350+425+428
!   - global reduction            --> end_method 1-3, lines 351-427,
!                                                           270-278
!   - different global            --> communication_method 1-7,
!     communication methods           lines 434-469, 518-809
 
 
!====================================================================== 

 
!CRAY: sn6715 hwwt3e 2.0.4.48 unicosmk CRAY T3E mpt.1.3.0.0.6 (514 nodes)
! 
!f90 -o heat-mpi1-big heat-mpi1-big.f 
!mpirun -np {size} ./heat-mpi1-big         (--> dedicated processors)
!
!size iter-  wall clock time communication part  abort criterion
!     ations    [seconds]    method [seconds]    meth. stride [seconds]
!
!   1  14320     10.01         0     0.000       1       1   0.4713E-01
!   1  14320     9.964         0     0.000       1      80   0.3292E-01
!   1  14406     10.04         0     0.000       1      98   0.3314E-01
!   1  14320     9.961         0     0.000       1     179   0.3271E-01
!
!   4  14320     4.945         1     2.536       1     179   0.1734E-01
!   4  14320     4.838         2     2.437       1     179   0.1505E-01
!   4  14320     5.005         3     2.574       1     179   0.1514E-01
!   4  14320     5.003         4     2.578       1     179   0.1495E-01
!   4  14320     4.305         5     1.878       1     179   0.1593E-01
!   4  14320     4.154         6     1.740       1     179   0.1423E-01
!   4  14320     2.858         7    0.4756       1     179   0.1378E-01
!   4  14320     5.239         5     2.023       1       1   0.7620    
!   4  14320     4.381         5     1.876       1      10   0.9422E-01
!   4  14320     4.289         5     1.824       1      20   0.5313E-01
!   4  14320     4.280         5     1.858       1      80   0.2137E-01
!   4  14406     4.247         5     1.813       1      98   0.1866E-01
!   4  14320     4.211         5     1.803       1     179   0.1423E-01
!   4  14321     4.186         5     1.787       2       0   0.7035E-02
!   4  14324     4.199         5     1.788       3       0   0.1827E-01
!   4  14320     3.709         7    0.6038       1       1   0.6805    
!   4  14320     2.858         7    0.4707       1      80   0.1955E-01
!   4  14406     2.875         7    0.4747       1      98   0.1757E-01
!   4  14320     2.850         7    0.4687       1     179   0.1368E-01
!   4  14320     2.862         7    0.4740       1      80   0.1937E-01
!   4  14406     2.873         7    0.4728       1      98   0.1758E-01
!   4  14320     2.850         7    0.4692       1     179   0.1373E-01
! 
!  16  14320     4.096         1     3.456       1     179   0.2059E-01
!  16  14320     3.993         2     3.353       1     179   0.2650E-01
!  16  14320     4.087         3     3.432       1     179   0.3633E-01
!  16  14320     4.008         4     3.357       1     179   0.3308E-01
!  16  14320     2.478         5     1.843       1     179   0.1842E-01
!  16  14320     2.548         6     1.915       1     179   0.1795E-01
!  16  14320     2.440         6     1.806       1     179   0.1819E-01
!  16  14320     1.450         7    0.8139       1     179   0.1868E-01
!  16  14320     3.619         5     1.952       1       1    1.022    
!  16  14320     2.709         5     1.950       1      10   0.1395    
!  16  14320     2.650         5     1.957       1      20   0.7442E-01
!  16  14320     2.600         5     1.957       1      80   0.2811E-01
!  16  14329     2.602         5     1.957       1      89   0.2809E-01
!  16  14320     2.572         5     1.938       1     179   0.1835E-01
!  16  13749     2.470         5     1.885       2       0   0.8576E-02
!  16  14330     2.622         5     1.991       3       0   0.1362E-01
!  16  14320     2.422         7    0.7134       1       1    1.071    
!  16  14320     1.456         7    0.8121       1      80   0.2739E-01
!  16  14406     1.461         7    0.8143       1      98   0.2476E-01
!  16  14320     1.444         7    0.8085       1     179   0.1852E-01
!  16  14320     1.457         7    0.8117       1      80   0.2762E-01
!  16  14406     1.458         7    0.8138       1      98   0.2514E-01
!  16  14320     1.443         7    0.8079       1     179   0.1856E-01
! 
!f90 -O3 -o heat-mpi1-big heat-mpi1-big.f 
!mpirun -np {size} ./heat-mpi1-big         (--> dedicated processors)
!
!size iter-  wall clock time communication part  abort criterion
!     ations    [seconds]    method [seconds]    meth. stride [seconds]
!
!   1  14320     9.964         0     0.000       1       1   0.4420E-01
!   1  14320     9.955         0     0.000       1      80   0.3276E-01
!   1  14406     10.01         0     0.000       1      98   0.3258E-01
!   1  14320     9.949         0     0.000       1     179   0.3282E-01
 
 
!====================================================================== 
 

!Hitachi: HI-UX/MPP hitachi 02-03 0 SR2201 (32 compute nodes)
! 
!mpif90 -o heat-mpi1-big heat-mpi1-big.f 
!mpirun -n {size} ./heat-mpi1-big         (--> dedicated processors)
! 
!size iter-  wall clock time communication part  abort criterion
!     ations    [seconds]    method [seconds]    meth. stride [seconds]
!
!   1  14320     26.48         0    0.0          1       1   0.7768D-01
!   1  14320     26.48         0    0.0          1      80   0.6874D-01
!   1  14406     26.63         0    0.0          1      98   0.6887D-01
!   1  14320     26.48         0    0.0          1     179   0.6789D-01
!
!   4  14320     10.77         1     3.827       1     179   0.6805D-01
!   4  14320     10.76         2     3.816       1     179   0.6786D-01
!   4  14320     10.60         3     3.640       1     179   0.7433D-01
!   4  14320     10.60         4     3.648       1     179   0.7391D-01
!   4  14320     11.40         5     4.465       1     179   0.6608D-01
!   4  14320     10.97         6     4.036       1     179   0.6573D-01
!   4  14320     11.52         7     4.578       1     179   0.6673D-01
!   4  14320     13.16         5     4.528       1       1    1.752    
!   4  14320     11.58         5     4.489       1      10   0.2267    
!   4  14320     11.49         5     4.476       1      20   0.1418    
!   4  14320     11.42         5     4.473       1      80   0.7830D-01
!   4  14406     11.48         5     4.496       1      98   0.7470D-01
!   4  14320     11.41         5     4.470       1     179   0.6594D-01
!   4  14321     11.40         5     4.482       2       0   0.5179D-01
!   4  14324     11.42         5     4.494       3       0   0.5899D-01
!   4  14320     13.24         7     4.537       1       1    1.827    
!   4  14320     11.52         7     4.573       1      80   0.7893D-01
!   4  14406     11.59         7     4.602       1      98   0.7541D-01
!   4  14320     11.51         7     4.575       1     179   0.6683D-01
!   4  14320     11.52         7     4.573       1      80   0.7906D-01
!   4  14406     11.59         7     4.602       1      98   0.7555D-01
!   4  14320     11.51         7     4.574       1     179   0.6644D-01
! 
!  16  14320     7.196         1     5.338       1     179   0.8017D-01
!  16  14320     7.360         2     5.498       1     179   0.8334D-01
!  16  14320     7.102         3     5.214       1     179   0.9867D-01
!  16  14320     7.230         4     5.343       1     179   0.1003    
!  16  14320     7.748         5     5.889       1     179   0.7961D-01
!  16  14320     7.126         6     5.272       1     179   0.7837D-01
!  16  14320     7.122         6     5.269       1     179   0.7831D-01
!  16  14320     17.40         7     15.51       1     179   0.7954D-01
!  16  14320     10.81         5     4.819       1       1    4.191    
!  16  14320     8.027         5     5.766       1      10   0.4789    
!  16  14320     7.872         5     5.823       1      20   0.2680    
!  16  14320     7.760         5     5.870       1      80   0.1093    
!  16  14329     7.756         5     5.871       1      89   0.1039    
!  16  14320     7.735         5     5.875       1     179   0.8000D-01
!  16  13749     7.393         5     5.680       2       0   0.4482D-01
!  16  14330     7.840         5     6.006       3       0   0.5258D-01
!  16  14320     20.49         7     15.26       1       1    3.405    
!  16  14320     17.37         7     15.46       1      80   0.1019    
!  16  14406     17.46         7     15.55       1      98   0.9569D-01
!  16  14320     17.35         7     15.46       1     179   0.7946D-01
!  16  14320     17.37         7     15.46       1      80   0.1030    
!  16  14406     17.46         7     15.55       1      98   0.9452D-01
!  16  14320     17.34         7     15.46       1     179   0.7895D-01

!f90 -W0,'opt(o(ss))' -I/usr/local/mpi/include heat-mpi1-big.f -o heat-mpi1-big  -L/usr/local/mpi/lib/hmpp2/cml -lfmpi -lmpi
!mpirun -n {size} ./heat-mpi1-big         (--> dedicated processors)
! 
!size iter-  wall clock time communication part  abort criterion
!     ations    [seconds]    method [seconds]    meth. stride [seconds]
!
!   1  14320     21.17         0    0.0          1       1   0.5596D-01
!   1  14320     21.17         0    0.0          1      80   0.5276D-01
!   1  14406     21.29         0    0.0          1      98   0.5291D-01
!   1  14320     21.17         0    0.0          1     179   0.5157D-01
 
 
!====================================================================== 


!NEC: SUPER-UX hwwsx4 8.1 Rev1 SX-4 (32 processors)
! 
!mpif90 -o heat-mpi1-big heat-mpi1-big.f 
!mpirun -np {size} ./heat-mpi1-big          (--> time sharing)
!
!size iter-  wall clock time communication part  abort criterion
!     ations    [seconds]    method [seconds]    meth. stride [seconds]
!
!   1  14320     2.254         0    0.0000D+00   1       1   0.7463D-02
!   1  14320     2.407         0    0.0000D+00   1      80   0.2452D-01
!   1  14406     3.248         0    0.0000D+00   1      98   0.1007D-01
!   1  14320     3.275         0    0.0000D+00   1     179   0.1284D-01
!
!   4  aborting in communication1 
!   4  aborting in communication2 
!   4  14320     5.797         3     4.515       1     179   0.7583D-01
!   4  14320     9.017         4     7.778       1     179   0.8718D-01
!   4  14320     8.053         5     6.862       1     179   0.8165D-01
!   4  14320     4.110         6     3.008       1     179   0.2188D-01
!   4  14320     7.352         7     6.111       1     179   0.3084D-01
!   4  14320     8.788         5     4.794       1       1    2.668    
!   4  14320     7.327         5     5.756       1      10   0.4210    
!   4  14320     10.63         5     9.075       1      20   0.3413    
!   4  14320     4.882         5     3.736       1      80   0.6535D-01
!   4  14406     4.497         5     3.377       1      98   0.3326D-01
!   4  14320     5.314         5     4.153       1     179   0.2782D-01
!   4  14321     5.331         5     4.243       2       0   0.2061D-01
!   4  14324     7.461         5     6.184       3       0   0.4413D-01
!   4  14320     9.091         7     3.690       1       1    4.237    
!   4  14320     5.914         7     4.681       1      80   0.1237    
!   4  14406     4.833         7     3.658       1      98   0.3768D-01
!   4  14320     7.123         7     5.976       1     179   0.5228D-01
!   4  14320     10.29         7     8.833       1      80   0.1718    
!   4  14406     7.971         7     6.813       1      98   0.3733D-01
!   4  14320     5.365         7     4.205       1     179   0.6619D-01
 
 
!====================================================================== 


!HP: HP-UX hp-v B.11.00 A 9000/800 75859 two-user license (8 processors)
!
!mpif90 -o heat-mpi1-big heat-mpi1-big.f
!mpirun -np {size} heat-mpi1-big           (--> time sharing)
! 
!size iter-  wall clock time communication part  abort criterion
!     ations    [seconds]    method [seconds]    meth. stride [seconds]
!
!   1  14320     154.2         0     0.000       1       1   0.1943    
!   1  14320     153.8         0     0.000       1      80   0.1912    
!
!   4  14320     40.93         1     2.829       1     179   0.1513    
!   4  14320     40.83         2     3.112       1     179   0.1480    
!   4  14320     40.81         3     2.928       1     179   0.1496    
!   4  14320     40.59         4     2.688       1     179   0.1584    
!   4  14320     40.67         5     2.920       1     179   0.1407    
!   4  14320     40.17         6     2.399       1     179   0.1441    
!   4  15036     42.37         7     2.423       1     179   0.1574    
!   4  14320     41.01         5     2.482       1       1   0.8361    
!   4  14320     40.99         5     2.938       1      10   0.2516    
!   4  14320     40.84         5     3.089       1      20   0.1919    
!   4  14320     41.42         5     3.248       1      80   0.1553    
!   4  14406     41.48         5     3.272       1      98   0.1508    
!   4  14320     40.59         5     2.786       1     179   0.1408    
!   4  14321     41.00         5     3.143       2       0   0.1335    
!   4  14324     40.38         5     2.732       3       0   0.1313    
!   4  14521     41.50         7     1.598       1       1    1.414    
!   4  14640     41.40         7     2.367       1      80   0.1618    
!   4  14602     40.83         7     2.191       1      98   0.1545    
!   4  15036     42.19         7     2.336       1     179   0.1549    
!   4  14602     41.47         7     2.837       1      98   0.1642    
!   4  15036     42.44         7     2.652       1     179   0.1629    
!
!   8  16468     23.83         1     2.853       1     179   0.2768    
!   8  13067     18.79         2     2.145       1     179   0.2052    
!   8  13067     19.14         3     2.579       1     179   0.1634    
!   8  13067     18.46         4     1.919       1     179   0.1861    
!   8  13067     18.88         5     2.211       1     179   0.2258    
!   8  13067     18.48         6     1.828       1     179   0.2099    
!   8  14141     20.41         7     2.421       1     179   0.2183    
!   8  12973     20.65         5     1.815       1       1    1.806    
!   8  12980     20.86         5     3.062       1      10   0.5012    
!   8  12980     20.20         5     2.550       1      20   0.2519    
!   8  13040     25.13         5     7.000       1      80   0.3170    
!   8  13034     18.64         5     2.051       1      98   0.1903    
!   8  13067     20.51         5     3.744       1     179   0.2561    
!   8  12789     18.24         5     2.262       2       0   0.1765    
!   8  12980     19.35         5     2.774       3       0   0.1791    
!   8  13102     19.63         7     1.726       1       1    1.307    
!   8  13520     19.34         7     2.127       1      80   0.2165    
!   8  13230     19.81         7     2.837       1      98   0.2056    
!   8  14141     20.30         7     2.307       1     179   0.2041    
! 
!mpif90 -O3 -o heat-mpi1-big heat-mpi1-big.f
!mpirun -np {size} heat-mpi1-big           (--> time sharing)
! 
!size iter-  wall clock time communication part  abort criterion
!     ations    [seconds]    method [seconds]    meth. stride [seconds]
!
!   1  14320     7.946         0     0.000       1       1   0.1201    
!   1  14320     7.934         0     0.000       1      80   0.1193    
!   1  14406     7.981         0     0.000       1      98   0.1201    
!   1  14320     7.918         0     0.000       1     179   0.1192    
!
!   4  14320     4.873         1     2.211       1     179   0.2076    
!   4  14320     5.059         2     2.341       1     179   0.2096    
!   4  14320     4.624         3     2.059       1     179   0.1613    
!   4  14320     4.394         4     1.775       1     179   0.1872    
!   4  14320     4.642         5     2.141       1     179   0.1520    
!   4  14320     4.003         6     1.343       1     179   0.2064    
!   4  15036     3.986         7     1.230       1     179   0.2124    
!   4  14320     5.074         5     2.172       1       1   0.5087    
!   4  14320     4.900         5     2.152       1      10   0.2088    
!   4  14320     4.648         5     2.116       1      20   0.1789    
!   4  14320     4.634         5     2.127       1      80   0.1567    
!   4  14406     4.712         5     2.175       1      98   0.1677    
!   4  14320     4.655         5     2.124       1     179   0.1678    
!   4  14321     4.722         5     2.111       2       0   0.1643    
!   4  14324     4.681         5     2.130       3       0   0.1572    
!   4  14521     4.500         7     1.248       1       1   0.7473    
!   4  14640     3.982         7     1.282       1      80   0.2118    
!   4  14602     3.875         7     1.194       1      98   0.2069    
!   4  15036     4.037         7     1.270       1     179   0.2069    
!   4  14640     4.094         7     1.337       1      80   0.2020    
!   4  14602     3.948         7     1.248       1      98   0.2036    
!   4  15036     4.056         7     1.302       1     179   0.2108    
 
 
!====================================================================== 


!SGI: IRIX64 vision 6.4 02121744 IP27  (14 R10000 CPUs, distributed)
! 
!f90 -o heat-mpi1-big heat-mpi1-big.f -lmpi
!mpirun -np 1 ./heat-mpi1-big             (--> time sharing)
!
!size iter-  wall clock time communication part  abort criterion
!     ations    [seconds]    method [seconds]    meth. stride [seconds]
!
!   1  14320     65.09         0     0.000       1       1   0.9623E-02
!   1  14320     65.07         0     0.000       1      80   0.8041E-02
! 
!f90 -O3 -o heat-mpi1-big heat-mpi1-big.f -lmpi
!mpirun -np {size} ./heat-mpi1-big        (--> time sharing)
!
!size iter-  wall clock time communication part  abort criterion
!     ations    [seconds]    method [seconds]    meth. stride [seconds]
!
!   1  14320     5.685         0     0.000       1       1   0.9184E-02
!   1  14320     5.678         0     0.000       1      80   0.8256E-02
!   1  14406     5.714         0     0.000       1      98   0.8275E-02
!   1  14320     5.681         0     0.000       1     179   0.1106E-01
!
!   4  14320     3.326         1     2.119       1     179   0.1339E-01
!   4  14320     3.243         2     2.056       1     179   0.1310E-01
!   4  14320     2.955         3     1.775       1     179   0.1345E-01
!   4  14320     2.894         4     1.708       1     179   0.1275E-01
!   4  14320     3.271         5     2.094       1     179   0.1227E-01
!   4  14320     2.811         6     1.642       1     179   0.1280E-01
!   4  14320     2.940         7     1.755       1     179   0.1296E-01
!   4  14320     3.999         5     2.187       1       1   0.6505    
!   4  14320     3.315         5     2.084       1      10   0.7362E-01
!   4  14320     3.288         5     2.082       1      20   0.4114E-01
!   4  14320     3.228         5     2.058       1      80   0.1642E-01
!   4  14406     3.301         5     2.112       1      98   0.1512E-01
!   4  14320     3.233         5     2.061       1     179   0.1191E-01
!   4  14321     3.231         5     2.082       2       0   0.9707E-02
!   4  14324     3.190         5     2.058       3       0   0.1644E-01
!   4  14320     3.956         7     1.796       1       1    1.020    
!   4  14320     2.936         7     1.778       1      80   0.2272E-01
!   4  14406     2.950         7     1.787       1      98   0.1996E-01
!   4  14320     2.930         7     1.778       1     179   0.1424E-01
!   4  14320     2.934         7     1.775       1      80   0.2260E-01
!   4  14406     2.953         7     1.788       1      98   0.2116E-01
!   4  14320     2.931         7     1.778       1     179   0.1448E-01
!
!   8 --> deadlock with all communication_methods !!!!!!!!!!!!!!! 
!
!   9  14320     3.500         1     2.951       1     179   0.1654E-01
!   9  14320     3.389         2     2.842       1     179   0.1532E-01
!   9  14320     3.026         3     2.470       1     179   0.2244E-01
!   9  14320     2.997         4     2.439       1     179   0.2102E-01
!   9  14320     3.016         5     2.466       1     179   0.1622E-01
!   9  14320     2.722         6     2.172       1     179   0.1569E-01
!   9  14320     4.268         7     3.692       1     179   0.2562E-01
!   9  14320     4.770         5     2.601       1       1    1.632    
!   9  14320     3.211         5     2.528       1      10   0.1537    
!   9  14320     3.116         5     2.508       1      20   0.8018E-01
!   9  14320     3.064         5     2.508       1      80   0.2573E-01
!   9  14406     3.120         5     2.564       1      98   0.2244E-01
!   9  14320     3.051         5     2.502       1     179   0.1646E-01
!   9  14057     3.027         5     2.503       2       0   0.9316E-02
!   9  14326     3.119         5     2.574       3       0   0.1452E-01
!   9  14320     6.319         7     3.566       1       1    2.207    
!   9  14320     4.331         7     3.740       1      80   0.4852E-01
!   9  14406     4.360         7     3.771       1      98   0.4070E-01
!   9  14320     4.336         7     3.766       1     179   0.2584E-01
!   9  14320     4.353         7     3.762       1      80   0.4824E-01
!   9  14406     4.406         7     3.818       1      98   0.4005E-01
!   9  14320     4.335         7     3.763       1     179   0.2583E-01
 
 
!====================================================================== 
 
!Fujitsu VPP700: UNIX_System_V vpc004 4.1 ES 2 F700 UXP/V (52 vector PEs, distr.)
!frt heat-mpi1-big.f -Wl,-P,-J -dn -L/usr/lang/mpi/lib -lmpi -lmp -lpx -lelf -lgen -I/usr/lang/mpi/include -o heat-mpi1-big
!echo './heat-mpi1-big' > heat-mpi1-big.sh
! 
!qsub -lM 100mb -lT 1000 -lP 4 heat-mpi1-big.sh
!
!size iter-  wall clock time communication part  abort criterion
!     ations    [seconds]    method [seconds]    meth. stride [seconds]
!
!   4  14320     5.833         1     3.632       1     179   0.5760e-01
!   4  14320     5.950         2     3.761       1     179   0.5752e-01
!   4  14320     5.761         3     3.588       1     179   0.5455e-01
!   4  14320     6.047         4     3.868       1     179   0.5444e-01
!   4  14320     6.170         5     3.995       1     179   0.5338e-01
!   4  14320     5.767         6     3.595       1     179   0.5049e-01
!   4  14320     6.131         7     3.952       1     179   0.5608e-01
!   4  14320     8.797         5     4.708       1       1    1.932    
!   4  14320     6.412         5     4.055       1      10   0.2343    
!   4  14320     6.289         5     4.029       1      20   0.1406    
!   4  14320     6.214         5     4.016       1      80   0.6800e-01
!   4  14406     6.250         5     4.057       1      98   0.6316e-01
!   4  14320     6.165         5     3.991       1     179   0.5382e-01
!   4  14321     6.159         5     3.986       2       0   0.4633e-01
!   4  14324     6.189         5     3.970       3       0   0.8019e-01
!   4  14320     8.403         7     4.242       1       1    2.017    
!   4  14320     6.166         7     3.966       1      80   0.7835e-01
!   4  14406     6.190         7     3.985       1      98   0.7138e-01
!   4  14320     6.127         7     3.950       1     179   0.5625e-01
!   4  14320     6.168         7     3.963       1      80   0.7985e-01
!   4  14406     6.226         7     4.023       1      98   0.7029e-01
!   4  14320     6.159         7     3.981       1     179   0.5583e-01
! 
! 
!qsub -lM 100mb -lT 1000 -lP 8 heat-mpi1-big.sh
!
!size iter-  wall clock time communication part  abort criterion
!     ations    [seconds]    method [seconds]    meth. stride [seconds]
!
!   8  13246     5.330         1     4.225       1     179   0.6158e-01
!   8  13246     5.343         2     4.242       1     179   0.5888e-01
!   8  13246     4.897         3     3.786       1     179   0.6707e-01
!   8  13246     5.077         4     3.958       1     179   0.6669e-01
!   8  13246     5.538         5     4.446       1     179   0.5525e-01
!   8  13246     5.107         6     4.019       1     179   0.5360e-01
!   8  13246     7.502         7     6.389       1     179   0.7562e-01
!   8  13156     8.241         5     4.676       1       1    2.512    
!   8  13160     5.737         5     4.399       1      10   0.3037    
!   8  13160     5.559         5     4.355       1      20   0.1714    
!   8  13200     5.491         5     4.381       1      80   0.7372e-01
!   8  13230     5.455         5     4.347       1      98   0.6802e-01
!   8  13246     5.471         5     4.377       1     179   0.5457e-01
!   8  12971     5.364         5     4.322       2       0   0.4193e-01
!   8  13163     5.460         5     4.352       3       0   0.7398e-01
!   8  13156     10.05         7     6.087       1       1    2.907    
!   8  13200     7.484         7     6.326       1      80   0.1243    
!   8  13230     7.520         7     6.375       1      98   0.1083    
!   8  13246     7.545         7     6.433       1     179   0.7724e-01
!   8  13200     7.520         7     6.366       1      80   0.1215    
!   8  13230     7.491         7     6.345       1      98   0.1081    
!   8  13246     7.548         7     6.434       1     179   0.7727e-01
! 
! 
!qsub -lM 100mb -lT 1000 -lP 16 heat-mpi1-big.sh
!
!size iter-  wall clock time communication part  abort criterion
!     ations    [seconds]    method [seconds]    meth. stride [seconds]
!
!  16  14320     7.104         1     5.930       1     179   0.7022e-01
!  16  14320     7.146         2     5.971       1     179   0.7044e-01
!  16  14320     6.462         3     5.271       1     179   0.8718e-01
!  16  14320     6.687         4     5.500       1     179   0.8890e-01
!  16  14320     7.208         5     6.038       1     179   0.6749e-01
!  16  14320     6.474         6     5.313       1     179   0.6364e-01
!  16  14320     14.73         7     13.51       1     179   0.1249    
!  16  14320     10.84         5     5.706       1       1    4.015    
!  16  14320     7.585         5     6.029       1      10   0.4547    
!  16  14320     7.449         5     6.081       1      20   0.2700    
!  16  14320     7.328         5     6.132       1      80   0.9851e-01
!  16  14406     7.289         5     6.098       1      98   0.8718e-01
!  16  14320     7.273         5     6.106       1     179   0.6710e-01
!  16  13749     6.952         5     5.877       2       0   0.4500e-01
!  16  14330     7.444         5     6.249       3       0   0.8170e-01
!  16  14320     18.38         7     12.14       1       1    5.128    
!  16  14320     14.78         7     13.46       1      80   0.2279    
!  16  14406     14.90         7     13.61       1      98   0.1928    
!  16  14320     14.81         7     13.60       1     179   0.1238    
!  16  14320     14.73         7     13.40       1      80   0.2256    
!  16  14406     14.86         7     13.56       1      98   0.1923    
!  16  14320     14.68         7     13.46       1     179   0.1225    
! 
!---------------------------------------------------------------------- 
!
!heat-mpi2-big.sh :
!#!/bin/sh 
!LD_LIBRARY_PATH=/usr/lang/mpi2/lib:$LD_LIBRARY_PATH
!export LD_LIBRARY_PATH
!./heat-mpi2-big 
! 
!/usr/lang/mpi2/bin/mpifrt -o heat-mpi2-big heat-mpi1-big.f
!
!qsub -lM 100mb -lT 1000 -lP  4 heat-mpi2-big.sh
!
!size iter-  wall clock time communication part  abort criterion
!     ations    [seconds]    method [seconds]    meth. stride [seconds]
! 
!   4  14320     5.720         1     3.456       1     179   0.6613e-01
!   4  14320     5.681         2     3.423       1     179   0.6613e-01 
!   aborted in routine communication3 
!
!qsub -lM 100mb -lT 1000 -lP  8 heat-mpi2-big.sh
!
!
!qsub -lM 100mb -lT 1000 -lP 16 heat-mpi2-big.sh 
!
 

 
 
!====================================================================== 
