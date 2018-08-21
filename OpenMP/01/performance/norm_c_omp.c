/*
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
*/

#include <stdio.h>
#include <math.h>
#include "hlrs_get_time.h"

#ifdef _OPENMP
#   include <omp.h>
#endif

double a[4097][4097];

double norm(int n){
  double mymax[4097];
  int  i,j;
  double result=0.0;
  
#pragma omp parallel private(i,mymax)
  {
#pragma omp parallel for
    for(j=0; j<n; j++){
      mymax[j]=0;
      for(i=0; i<n ; i++)
	mymax[j] = mymax[j] + (a[j][i] >0.0 ? (a[j][i]) : (-a[j][i]) );
      /*	 mymax[j] = mymax[j] + fabs(a[j][i]);*/
    }

    result=0;
#pragma omp parallel for
    for(j=0;j<n; j++){
#pragma omp critical
      if (mymax[j]>result)
	result=mymax[j];
    }

  }/* end parallel */
  return result;
}
    

double norm_triangular(int n){
  double mymax;
  int  i,j;
  double result=0.0;
  
#pragma omp parallel for private(i,mymax) 
  for(j=0; j<n; j++){
    mymax=0;
    for(i=0; i<=j ; i++)
      mymax = mymax + (a[j][i] >0.0 ? (a[j][i]) : (-a[j][i]) );
    /*	 mymax = mymax + fabs(a[j][i]);*/
#pragma omp critical	  
    if ( mymax > result)
      result=mymax;
  }
  return result;
}
    

int main(int argc, char** argv){

  double reals,users;
  double reale,usere;
  double totals,totale;
  double s;


  int n,i,j,steps;
  int cpus, threads;


  hlrs_get_time(&totals,&users);
  cpus=0;
  threads=0;
#ifdef _OPENMP
#pragma omp parallel
  {
#pragma omp master
    {
      cpus = omp_get_num_procs();
      threads = omp_get_num_threads();
    }
  }
#endif


  printf("# OpenMP Version %d threads on %d CPUs \n",threads,cpus);
  printf("# size   real  repetitions Mops/s\n");

  for(i=0; i < 4097; i++)
    for(j=0; j < 4097 ; j++)
      a[i][j] = 1.0;

  n=2;
  for(steps=0;steps<12;steps++){
    i = 0;
    norm(n);
    hlrs_get_time(&users,&reals);
    do {
      i = i+1;
      s=norm(n);
      hlrs_get_time(&usere,&reale);
    } while ( usere-users < 0.2 && i<100000000);
    printf("%d %g %d %g\n", n,(reale-reals)/(double)(i),i,(double)(i)*n*n/(reale-reals));
    if (s != n) 
      printf("# wrong result!  %d != % g\n", n, s);
    n=n*2;

  }

  printf("#  Results of Triangular Matrix:\n");

  n=2;
  for(steps=0;steps<12;steps++){
    i = 0;
    norm_triangular(n);
    hlrs_get_time(&users,&reals);
    do {
      i = i+1;
      s=norm_triangular(n);
      hlrs_get_time(&usere,&reale);
    } while ( usere-users < 0.2 && i<100000000);
    printf("#%d %g %d %g\n", n,(reale-reals)/(double)(i),i,(double)(i)*n*n/(reale-reals));
    if (s != n) 
      printf("# wrong result!  %d != % g\n", n, s);
    n=n*2;

  }
  hlrs_get_time(&totale,&usere);
  printf("# total time = %g s\n",totale-totals);
  return 0;
}


