// FILE:
// heat-lame.c
//
// PORT:
// Stefan Kombrink
//
// DESCRIPTION:
// parallelized version of heat.c, a port from fortran heat-mpi-big0.f
// for details see heat-mpi-big0.f / heat-mpi-big0.h
//
// LAST MODIFIED:
// 18.01.2001
//
//




// for parallelization
#include <mpi.h>

// for text io
#include <stdio.h>
#include <time.h>

// basic maths
#include <math.h>

// for memory allocation
#include <stdlib.h>


// definitions for heat example
// DEFINES
#define FALSE 0
#define TRUE (!FALSE)
 
#define I     0 // for accessing the dimension of an array
#define K     1

// MACROS
#define minev(a, b) ( (a) < (b) ? (a) : (b) )
#define maxev(a, b) ( (a) > (b) ? (a) : (b) )

// for addressing the phi/phin array in the same way c addresses
// an equivalent two-dimensional array
#define addr(i, k) ( k-s[K] + (i-s[I])*(e[K]-s[K]+1) )
#define addrn(i, k) ( k-s[K]-b1 + (i-s[I]-b1)*(e[K]-s[K]-2*b1+1) )


// CONSTANTS
int start[2] = { 0,  0};
int   max[2] = {80, 80};
const int    itmax= 20000;
const int       b1=     1;

const double eps = 1.0e-08;


// GLOBALS
int upper, lower, right, left;
int dims[2], coords[2], s[2], e[2];
int mySize, myRank;
int inner[2], outer[2];

MPI_Comm myComm;



// PROCEDURES for heat example



// DESCRIPTION:
// distributes given number of nodes
// as clever as possible on each dimension
// and dumps the results

void DistributeNodes(void)
{
  // do NOT use more nodes than cells exists
  if ( (start[I]-max[I]+1) * (start[K]-max[K]+1) < mySize )
  {
    if (myRank == 0)
      printf("\nERROR:do NOT use more nodes than cells exists!\n");
    MPI_Finalize();
    exit(0);
  }

  dims[I] = dims[K] = 0;
  MPI_Dims_create(mySize, 2, dims);

  if (myRank == 0)
  {
    printf("\nusing %d nodes for calculation...\n", mySize);
    printf("cutting field [%d, %d] in %dx%d slices...\n",
           max[I], max[K], dims[I], dims[K]);
  }
}

// DESCRIPTION:
// calculates the boundaries of each 
// single area calculated by own process
// and dumps the results

void CalcBoundaries(void)
{
  int c;

  int inner1, n;

  for (c=0; c<2; c++)
  {
    inner1 = ((1 + max[c] - start[c]) - 2 * b1 - 1) / dims[c] + 1;
    n      =  (1 + max[c] - start[c]) - 2 * b1 - dims[c] * (inner1 - 1);

    if (coords[c] < n)
    {
      inner[c] = inner1;
      s[c]     = start[c] + coords[c] * inner[c];
    } else
    {
      inner[c] = inner1 - 1;
      s[c]     = start[c] + n * inner1 + (coords[c] - n) * inner[c];
    }

    outer[c] = inner[c] + 2 * b1;
    e[c]     = s[c] + outer[c] - 1;
  }

  // wait for ALL processes (should be printed out synced!)
  MPI_Barrier(myComm);

  if (myRank == 0) printf("\n synced with MPI_Barrier...\n");
  printf("\nnode %d has boundaries [(%d, %d), (%d, %d)]\n",
         myRank, s[I], s[K], e[I], e[K]);

  MPI_Barrier(myComm);
  return;
}



void Algorithm(int stride, int prt);


// main procedure
int main(argc, argv)

int argc;
char **argv;

{
  int period[2];


  // MPI initialization
  MPI_Init(&argc, &argv);

  // get information about own process and number of nodes
  MPI_Comm_size(MPI_COMM_WORLD, &mySize);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  // domain decomposition
  DistributeNodes();


  // initializing variables
  period[I] = FALSE;
  period[K] = FALSE;

  // creating a new MPI world!
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, period, TRUE, &myComm);

  // get new rank & coords
  MPI_Comm_rank(myComm, &myRank);
  MPI_Cart_coords(myComm, myRank, 2, coords);


  // computing s, e, inner, outer for both dimensions
  // calculating boundaries
  CalcBoundaries();


  // CALL HERE THE MAIN ALGORITHM!!!
  Algorithm(1, FALSE);


  MPI_Barrier(myComm);

  // wait for all processes to terminate
  MPI_Finalize();

  return 0;
}


void heatpr(double *);


void Algorithm(int stride, int prt)
{

  int i, k;
  int it;
 
  double *phi, *phin;
  double dx, dy, dx2, dy2, dx2i, dy2i;
  double dt, dphi, dphimax;

  double dphiMaxPartial;

  double startTime, endTime;
  double commTime, criterionTime;

  MPI_Status statuses[4];
  MPI_Request req[4];     
  MPI_Datatype horizontalBorder, verticalBorder;


  // allocate memory for dynamic arrays
  phi  = (double *)malloc( sizeof(double)*(e[K] - s[K] + 1)*(e[I] - s[I] + 1) );
  phin = (double *)malloc( sizeof(double)*(e[K] - s[K] - 2*b1 + 1)*(e[I] - s[I] - 2*b1 + 1));


  MPI_Cart_shift(myComm, I, TRUE, &left , &right);
  MPI_Cart_shift(myComm, K, TRUE, &lower, &upper);


  // create a MPI vector....

  MPI_Type_vector(b1, inner[K], outer[K], MPI_DOUBLE, &verticalBorder);
  MPI_Type_commit(&verticalBorder);
  MPI_Type_vector(inner[I], b1, outer[K], MPI_DOUBLE, &horizontalBorder);
  MPI_Type_commit(&horizontalBorder);

  // init...

  dx   = 1.0f / max[I];
  dy   = 1.0f / max[K];
  dx2  = dx * dx;     
  dy2  = dy * dy;
  dx2i = 1.0f / dx2; 
  dy2i = 1.0f / dy2;

  dt = minev(dx2, dy2) / 4.0;

  // start values 0.0

  for (k=s[K]; k<=minev(e[K], max[K]-1); k++)
    for (i=maxev(1, s[I]); i<=minev(e[I], max[I]-1); i++)
      phi[addr(i, k)] = 0.0;
  

  // start values 1.0

  if (e[K] == max[K])
    for (i=s[I]; i<=e[I]; i++)
      phi[addr(i, max[K])] = 1.0;


  // start values dx
  // calculates interpolated values for 
  // border regions

  if (s[I] == 0)
    for (k=s[K]; k<=e[K]; k++)
      phi[addr(0, k)] = (k-start[K])*dy;

  if (e[I] == max[I])
    for (k=s[K]; k<=e[K]; k++)
      phi[addr(max[I], k)] = (k-start[K])*dy;

  // print variables
  if (myRank == 0)
  {
    printf("\nheat conduction 2d\ndx=%f dy=%f dt=%f eps=%f\n",dx, dy, dt, eps);
  }


  // save current time
  startTime = MPI_Wtime();
  commTime  = 0;
  criterionTime = 0;

  // iteration

  for (it=1; it<=itmax; it++)
  {
    dphimax = 0.0;

    for (k=s[K]+b1; k<=e[K]-b1; k++)
      for (i=s[I]+b1; i<=e[I]-b1; i++)
      {
        dphi = (phi[addr(i+b1,k)]+phi[addr(i-b1,k)]-2*phi[addr(i,k)])*dx2i
             + (phi[addr(i,k+b1)]+phi[addr(i,k-b1)]-2*phi[addr(i,k)])*dy2i;
        dphi = dphi * dt;
        dphimax = maxev(dphimax, dphi);
        phin[addrn(i, k)] = phi[addr(i, k)] + dphi;
      }


    // save values
    for (k=s[K]+b1; k<=e[K]-b1; k++)
      for (i=s[I]+b1; i<=e[I]-b1; i++)
        phi[addr(i, k)] = phin[addrn(i, k)];


    // for optimization: allreduce only each stride's loop

    criterionTime = criterionTime - MPI_Wtime();

    if ((it % stride) == 0)
    {
      if (mySize > 1)
      {
        dphiMaxPartial = dphimax;
        MPI_Allreduce(&dphiMaxPartial, &dphimax, 1, MPI_DOUBLE, MPI_MAX, myComm);
      }

      // abort criterion
      if (dphimax < eps)
        goto endOfLoop;
    }

    criterionTime = criterionTime + MPI_Wtime();
    commTime = commTime - MPI_Wtime();


    // send and receive to/from upper/lower

    // receive horizontal borders (upper/lower)
    MPI_Irecv(&phi[addr(s[I]+b1, s[K])], 1, horizontalBorder,
              lower, MPI_ANY_TAG, myComm, req);
    MPI_Irecv(&phi[addr(s[I]+b1, e[K])], 1, horizontalBorder,
              upper, MPI_ANY_TAG, myComm, req + 1);

    // send horizontal borders (upper/lower)
    MPI_Isend(&phi[addr(s[I]+b1, e[K]-b1)], 1, horizontalBorder,
              upper, 0, myComm, req + 2);
    MPI_Isend(&phi[addr(s[I]+b1, s[K]+b1)], 1, horizontalBorder,
              lower, 0, myComm, req + 3);
    MPI_Waitall(4, req, statuses);

    // send and receive to / from left / right

    // receive vertical borders (left/right)
    MPI_Irecv(&phi[addr(s[I], s[K]+b1)], 1, verticalBorder,
              left, MPI_ANY_TAG, myComm, req);
    MPI_Irecv(&phi[addr(e[I], s[K]+b1)], 1, verticalBorder,
              right, MPI_ANY_TAG, myComm, req + 1);

    // send vertical borders (left/right)
    MPI_Isend(&phi[addr(e[I]-b1, s[K]+b1)], 1, verticalBorder,
              right, 0, myComm, req + 2);
    MPI_Isend(&phi[addr(s[I]+b1, s[K]+b1)], 1, verticalBorder,
              left, 0, myComm, req + 3);
    MPI_Waitall(4, req, statuses);

    commTime = commTime + MPI_Wtime();    
  }    

endOfLoop:

  // because it's ONE iteration too much!
  it--;

  criterionTime = criterionTime + MPI_Wtime();
  endTime = MPI_Wtime();

  // print array?
  if (prt) heatpr(phi);

  if (myRank==0)
  {

    printf("\nsize iter-   wall clock time    comm   part     abort  criterion");
    printf("\n     ations    [seconds]      method [seconds]   meth.  stride [seconds]");
    printf("\n%4d %6d  %3.12f    %1d   %3.9f    %d     %d  %3.12f\n",
           mySize, it, endTime - startTime, 1, commTime, 1, stride, criterionTime);

    printf("\ntime: %4.4f ms\n",(endTime-startTime)*1000.0);
  }

  return;
}




void heatpr(double *phi)
{
  int ci, ck, i, k;
  int firsti, lasti;
  int firstk, lastk;

  printf("\nprinting results of node %d\n", myRank);

  for (ci=0; ci<=dims[I]; ci++)
  {
    firsti = s[I] + b1;
    lasti  = e[I] - b1;

    if (ci == 0) firsti = firsti - b1;
    if (ci == dims[I]-1) lasti = lasti + b1;

    for (ck=0; ck<=dims[K]; ck++)
    {
      if ( (ci == coords[I]) && (ck == coords[K]) )
      {
        if (ck == 0)
          printf("\nprinting the %dth horizontal block...\n", ci);
      
        firstk = s[K] + b1;
        lastk  = e[K] - b1;

        if (ck == 0) firstk = firstk - b1;
        if (ck == dims[K]-1) lastk = lastk + b1;
       
        for (k=firstk; k<=lastk; k++)
	{
          printf("[%2d:%2d] ", myRank, k);
          for (i=firsti; i<=lasti; i++)
            printf(" %1.2f", phi[addr(i, k)]);
          printf("\n");
	}
     
      }
      MPI_Barrier(myComm);
    } // end for ck
  } // end for ci
  return;
}









