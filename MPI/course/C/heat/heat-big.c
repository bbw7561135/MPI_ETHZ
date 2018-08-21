#include <stdio.h>
#include <time.h>

#define min(A,B) ((A) < (B) ? (A) : (B))
#define max(A,B) ((A) > (B) ? (A) : (B))

const int imax = 80;
const int kmax = 80;
const int itmax = 20000;

int main()
{
  double eps = 1.0e-08;
  double phi[imax+1][kmax+1], phin[imax][kmax];
  double dx,dy,dx2,dy2,dx2i,dy2i,dt,dphi,dphimax;
  int i,k,it;

  clock_t t;

  dx=1.0/kmax;
  dy=1.0/imax;
  dx2=dx*dx;
  dy2=dy*dy;
  dx2i=1.0/dx2;
  dy2i=1.0/dy2;
  dt=min(dx2,dy2)/4.0;
/* start values 0.d0 */
  for (k=0;k<kmax;k++)
  {
    for (i=1;i<imax;i++)
    {
      phi[i][k]=0.0;
    }
  }
/* start values 1.d0 */
  for (i=0;i<=imax;i++)
  {
    phi[i][kmax]=1.0;
  }
/* start values dx */
  phi[0][0]=0.0;
  phi[imax][0]=0.0;
  for (k=1;k<kmax;k++)
  {
    phi[0][k]=phi[0][k-1]+dx;
    phi[imax][k]=phi[imax][k-1]+dx;
  }

  printf("\nHeat Conduction 2d\n");
  printf("\ndx = %12.4g, dy = %12.4g, dt = %12.4g, eps = %12.4g\n",
         dx,dy,dt,eps);

  t=clock();

/* iteration */
  for (it=1;it<=itmax;it++)
  {
    dphimax=0.;
    for (k=1;k<kmax;k++)
    {
      for (i=1;i<imax;i++)
      {
        dphi=(phi[i+1][k]+phi[i-1][k]-2.*phi[i][k])*dy2i
            +(phi[i][k+1]+phi[i][k-1]-2.*phi[i][k])*dx2i;
        dphi=dphi*dt;
        dphimax=max(dphimax,dphi);
        phin[i][k]=phi[i][k]+dphi;
      }
    }
/* save values */
    for (k=1;k<kmax;k++)
    {
      for (i=1;i<imax;i++)
      {
	phi[i][k]=phin[i][k];
      }
    }
    if(dphimax<eps) break;
  }

  t=clock();
  printf("\n%d iterations\n",it);
  printf("\nCPU time = %#12.4g sec\n",t/1000000.0);

}
