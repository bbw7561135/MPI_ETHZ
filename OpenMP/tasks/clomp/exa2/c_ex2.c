/* compute the potential energy of a collection of */
/* particles interacting via pairwise potential    */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <omp.h>

int n;
int niter;
double *xdim, *ydim;
 
double computePot();
void power(double *x, double*y, int i, int j);
void allocPositions();
void initPositions();
void updatePositions();
 
int main() {
   int i;
   double pot;

   n = 2000;
   niter = 10;

   allocPositions();
   initPositions();
   updatePositions();
 
   for( i=0; i<niter; i++ ) {
      pot = computePot();
      printf("Potential: %f\n", pot);
      updatePositions();
   }
}


void allocPositions()
{
   double * r_tmp2;
   int i;

   r_tmp2 = (double *)malloc(sizeof(double)*2*n);
   if (!r_tmp2)
   {
      printf("Error allocating memory\n");
      exit(-1);
   }

   xdim = r_tmp2;
   r_tmp2 += n;
   ydim = r_tmp2;
}
 
void initPositions() {
   int i, j;
 
   for( j=0; j<n; j++ )
   {
      xdim[j] = 0.5 + ( (double) rand() / (double) RAND_MAX );
      ydim[j] = 0.5 + ( (double) rand() / (double) RAND_MAX );
   }
}
 
 
void updatePositions() {
   int i, j;
 
   for( j=0; j<n; j++ )
   {
      xdim[j] -= 0.5 + ( (double) rand() / (double) RAND_MAX );
      ydim[j] -= 0.5 + ( (double) rand() / (double) RAND_MAX );
   }
}
 
 
double computePot() {
   int i, j;
   double pot = 0.0;
   double distx, disty, dist;
 
   #pragma omp parallel private(distx, disty, dist, i, j) 
   {
      #pragma omp for schedule(static) reduction(+: pot) 
      for( i=0; i<n; i++ ) {
	 for( j=0; j<i; j++ ) {
	    power(&distx, &disty, i, j);
	    dist = sqrt( distx + disty );
	    pot += 1.0 / dist;
	 }
      }
   }
   return pot;
}

void power(double *x, double*y, int i, int j)
{
    *x = pow((xdim[j] - xdim[i]), 2);
    *y = pow((ydim[j] - ydim[i]), 2);
}
