This files in the subdirectories have been written 
as sample solutions to an exercise in a course given 
at the High Performance Computing Centre Stuttgart (HLRS).
The examples in Ch2-Ch12 are based on the examples in the
MPI course of the Edinburgh Parallel Computing Centre (EPCC).
The examples in the other directories are added by HLRS.
It is made freely available with the understanding that
every copy of the files must include their header and that
HLRS and EPCC take no responsibility for the use of the
enclosed teaching material.

Authors: Joel Malard, Alan Simpson,            (EPCC)
         Rolf Rabenseifner, Traugott Streicher (HLRS)

Contact: rabenseifner@hlrs.de

Naming of the files and subdirectories

- 11: Examples should run with MPI-1.1.
      In Fortran, "INCLUDE 'mpif.h'" is used.
      The example may contain MPI routines that are 
      deprecated or removed in later versions of MPI.

- 20: Examples should run with MPI-2.0 and later.
      In Fortran, "USE mpi" is used, and
      may contain calls to MPI that were
      introduced in MPI-3.0

- 30: Examples require MPI-3.0 and later.
      In Fortran, "USE mpi_f08" is used.

- No such numbers (in C only):
      Valid since MPI-1.1. 
      Should not contain calls to MPI routines
      that are deprecated or removed in a later
      version of MPI.

Ch14 is a link to an old-style Fortran directory F/mpi_io.
This mpi_io part is not yet updated to MPI-3.0.
