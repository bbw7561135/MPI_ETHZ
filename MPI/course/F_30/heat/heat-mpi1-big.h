        include 'mpif.h' 

        parameter (imax=80,kmax=80)
        parameter (istart=0,kstart=0) 
        parameter (itmax=20000)
        double precision eps
        parameter (eps=1.d-08)

        integer        neighbor(4)
        integer        upper, lower, right, left, rank
        equivalence   (neighbor(1), upper) 
        common /ranks/ upper, lower, right, left, rank
c                    = ---------neighbor--------  rank 

        integer        idim, kdim, icoord, kcoord
        common /coord/ idim, kdim, icoord, kcoord

        integer        sts_array(MPI_STATUS_SIZE,8) 
        integer        sts_upper(MPI_STATUS_SIZE) 
        integer        sts_lower(MPI_STATUS_SIZE) 
        integer        sts_right(MPI_STATUS_SIZE) 
        integer        sts_left (MPI_STATUS_SIZE) 
        integer        sts_recvs(MPI_STATUS_SIZE,4) 
        integer        sts_sends(MPI_STATUS_SIZE,4) 
        equivalence   (sts_array(1,1),sts_upper(1)) 
        equivalence   (sts_array(1,1),sts_recvs(1,1)) 
        common /sts/   sts_upper,sts_lower,sts_right,sts_left,sts_sends 
c                   =  --------------sts_recvs--------------- sts_sends
c                   =  ------------------sts_array--------------------- 

        integer        size, comm, horizontal_border, vertical_border
        common /misc/  size, comm, horizontal_border, vertical_border
 
        integer        b1
        parameter     (b1=1)
 
        integer        kinner, kouter, ks, ke, iinner, iouter, is, ie 
        common /block/ kinner, kouter, ks, ke, iinner, iouter, is, ie 

        integer        lower_left_tag, upper_right_tag
        common /tags/  lower_left_tag, upper_right_tag
 
        integer        ierror 
