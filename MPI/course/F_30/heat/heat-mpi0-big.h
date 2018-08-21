        include 'mpif.h' 

        parameter (imax=80,kmax=80)
        parameter (istart=0,kstart=0) 
        parameter (itmax=20000)
        double precision eps
        parameter (eps=1.d-08)

        integer        upper, lower, right, left, rank
        common /ranks/ upper, lower, right, left, rank

        integer        idim, kdim, icoord, kcoord
        common /coord/ idim, kdim, icoord, kcoord

        integer        size, comm
        common /misc/  size, comm
 
        integer        b1
        parameter     (b1=1)
 
        integer        kinner, kouter, ks, ke, iinner, iouter, is, ie 
        common /block/ kinner, kouter, ks, ke, iinner, iouter, is, ie 

        integer        ierror 
