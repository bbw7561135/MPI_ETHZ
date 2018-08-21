/* Program usage:  mpirun -np <procs> ./heat_petsc [-help] [all PETSc options] */ 
                                             /* valid for PETSC library  2.2.0 */
static char help[] = "Compute steady temperature distribution for given \
temperatures on a boundary in parallel with KSP.\n\
Input parameters include:\n\
  -random_exact_sol : use a random exact solution vector\n\
  -view_exact_sol   : write exact solution vector to stdout\n\
  -view_sol_serial  : write solution grid to stdout (1 item per line)\n\
  -view_sol         : write solution grid to stdout (as matrix)\n\
  -view_sol_x  -draw_pause 3 : view solution  x  on a X window\n\
  -view_mat_x  -draw_pause 3 : view the matrix A on a X window\n\
  -m <mesh_x>       : number of mesh points in x-direction\n\
  -n <mesh_n>       : number of mesh points in y-direction\n\n";
/*T
   Concepts: KSP^basic parallel example;
   Concepts: KSP^Laplacian, 2d
   Concepts: Laplacian, 2d
   Processors: n
T*/
/* 
  Include "petscsles.h" so that we can use KSP solvers.  Note that this file
  automatically includes:
     petsc.h       - base PETSc routines   petscvec.h - vectors
     petscsys.h    - system routines       petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners
*/
#include "petscksp.h"
#include "petscda.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  Vec         x,b,u;  /* approx solution, RHS, exact solution */
  Mat         A;        /* linear system matrix */
  KSP         ksp;      /* linear solver context */
  PetscRandom rctx;     /* random number generator context */
  PetscReal   norm;     /* norm of solution error */
  int         i,j,I,J,Istart,Iend,m = 4,n = 4,its;
  PetscTruth  flg, flg1;
  PetscScalar v,h,one = 1.0,neg_one = -1.0;
  PC pc;  PCType pctype;  KSPType ksptype;
    
  PetscInitialize(&argc,&args,(char *)0,help);
  PetscOptionsGetInt(PETSC_NULL,"-m",&m,PETSC_NULL);
  PetscOptionsGetInt(PETSC_NULL,"-n",&n,PETSC_NULL);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
         Compute the matrix and right-hand-side vector that define
         the linear system, Ax = b.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /* 
     Create parallel matrix, specifying only its global dimensions.
     When using MatCreate(), the matrix format can be specified at
     runtime. Also, the parallel partitioning of the matrix is
     determined by PETSc at runtime.

     Performance tuning note:  For problems of substantial size,
     preallocation of matrix memory is crucial for attaining good 
     performance.  Since preallocation is not possible via the generic
     matrix creation routine MatCreate(), we recommend for practical 
     problems instead to use the creation routine for a particular matrix
     format, e.g.,
         MatCreateMPIAIJ() - parallel AIJ (compressed sparse row)
         MatCreateMPIBAIJ() - parallel block AIJ
     See the matrix chapter of the users manual for details.
  */
  MatCreate(PETSC_COMM_WORLD,&A); MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m*n,m*n); /*2.3.0*/
/*MatCreate(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,m*n,m*n,&A);   2.2.1 */
  MatSetFromOptions(A);
  /* 
     Currently, all PETSc parallel matrix formats are partitioned by
     contiguous chunks of rows across the processors.  Determine which
     rows of the matrix are locally owned. 
  */
  MatGetOwnershipRange(A,&Istart,&Iend);

  /* 
     Set matrix elements for the 2-D, five-point stencil in parallel.
      - Each processor needs to insert only elements that it owns
        locally (but any non-local elements will be sent to the
        appropriate processor during matrix assembly). 
      - Always specify global rows and columns of matrix entries.

     Note: this uses the less common natural ordering that orders first
     all the unknowns for x = h then for x = 2h etc; Hence you see J = I +- n
     instead of J = I +- m as you might expect. The more standard ordering
     would first do all variables for y = h, then y = 2h etc.

   */
  for (I=Istart; I<Iend; I++) { 
    v = -1.0; i = I/n; j = I - i*n;  
    if (i>0)   {J = I - n; MatSetValues(A,1,&I,1,&J,&v,INSERT_VALUES);}
    if (i<m-1) {J = I + n; MatSetValues(A,1,&I,1,&J,&v,INSERT_VALUES);}
    if (j>0)   {J = I - 1; MatSetValues(A,1,&I,1,&J,&v,INSERT_VALUES);}
    if (j<n-1) {J = I + 1; MatSetValues(A,1,&I,1,&J,&v,INSERT_VALUES);}
    v = 4.0; MatSetValues(A,1,&I,1,&I,&v,INSERT_VALUES);
  }

  /* 
     Assemble matrix, using the 2-step process:
       MatAssemblyBegin(), MatAssemblyEnd()
     Computations can be done while messages are in transition
     by placing code between these two statements.
  */
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

  /* 
     Create parallel vectors.
      - We form 1 vector from scratch and then duplicate as needed.
      - When using VecCreate(), VecSetSizes and VecSetFromOptions()
        in this example, we specify only the
        vector's global dimension; the parallel partitioning is determined
        at runtime. 
      - When solving a linear system, the vectors and matrices MUST
        be partitioned accordingly.  PETSc automatically generates
        appropriately partitioned matrices and vectors when MatCreate()
        and VecCreate() are used with the same communicator.  
      - The user can alternatively specify the local vector and matrix
        dimensions when more sophisticated partitioning is needed
        (replacing the PETSC_DECIDE argument in the VecSetSizes() statement
        below).
  */
  VecCreate(PETSC_COMM_WORLD,&u);
  VecSetSizes(u,PETSC_DECIDE,m*n);
  VecSetFromOptions(u);
  VecDuplicate(u,&b); 
  VecDuplicate(b,&x);

  /* 
     If 'random_exact_sol' is False set rhs vector and exact solution
     vector for boundary conditions:
       0.......1
       .       .
       .       . j = -1 .. n
       .       .
       0.......1
      i = -1 .. m
     To be computed: i=0..m-1,j=0..n-1
     If 'random_exact_sol' is True set random exact solution; 
     then compute right-hand-side vector.
  */
  PetscOptionsHasName(PETSC_NULL,"-random_exact_sol",&flg);
  if ( !flg) { 
    VecGetOwnershipRange(b,&Istart,&Iend); 
    h = 1.0/(m+1);
    for (I=Istart; I<Iend; I++) {
      v = 0;  i = I/n;  j = I - i*n;
      if (i==0)   v = v + /*u(-1,j):*/  h * 0;     /* = 0 */
      if (i==m-1) v = v + /*u(m, j):*/  h * (m+1); /* = 1 */
      if (j==0)   v = v + /*u(i,-1):*/  h * (i+1);
      if (j==n-1) v = v + /*u(i, n):*/  h * (i+1);
      if (v != 0) ; VecSetValues(b,1,&I,&v,INSERT_VALUES);
      v = /*u(i, j):*/ h * (i+1); VecSetValues(u,1,&I,&v,INSERT_VALUES);
    }
    VecAssemblyBegin(b); VecAssemblyBegin(u);
    VecAssemblyEnd(b); VecAssemblyEnd(u);
  } else { 
    PetscRandomCreate(PETSC_COMM_WORLD,&rctx); 
    VecSetRandom(u,rctx);  PetscRandomDestroy(rctx);
    MatMult(A,u,b); 
  } 

  /*
     View the exact solution vector if desired
  */
  PetscOptionsHasName(PETSC_NULL,"-view_exact_sol",&flg);
  if (flg) {VecView(u,PETSC_VIEWER_STDOUT_WORLD);}

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                Create the linear solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* 
     Create linear solver context
  */
  KSPCreate(PETSC_COMM_WORLD,&ksp);

  /* 
     Set operators. Here the matrix that defines the linear system
     also serves as the preconditioning matrix.
  */
  KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);

  /* 
     Set linear solver defaults for this problem (optional).
     - By extracting the PC context from the KSP context,
       we can then directly call any KSP and PC routines to set
       various options.
     - The following two statements are optional; all of these
       parameters could alternatively be specified at runtime via
       KSPSetFromOptions().  All of these defaults can be
       overridden at runtime, as indicated below.
  */

  KSPSetTolerances(ksp,1.e-2/((m+1)*(n+1)),1.e-50,PETSC_DEFAULT,PETSC_DEFAULT);

  /* 
    Set runtime options, e.g.,
        -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    These options will override those specified above as long as
    KSPSetFromOptions() is called _after_ any other customization
    routines.
  */
  KSPSetFromOptions(ksp);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                      Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* KSPSetRhs(ksp,b); KSPSetSolution(ksp,x); KSPSolve(ksp); -* PETSc 2.2.0 */
  KSPSolve(ksp,b,x);                                         /* PETSc 2.2.1 */

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                      Print used solver (KSP and PC)
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 
  KSPGetType(ksp, &ksptype);
  KSPGetPC(ksp, &pc);
  PCGetType(pc, &pctype);
  PetscPrintf( PETSC_COMM_WORLD, "Used  -ksp_type %s   -pc_type %s \n", ksptype, pctype );

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                      Check solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* 
     Draw solution grid 
  */

  PetscOptionsHasName(PETSC_NULL,"-view_sol_serial",&flg);
  if (flg) {VecView(x,PETSC_VIEWER_STDOUT_WORLD);}
 
  PetscOptionsHasName(PETSC_NULL,"-view_sol",&flg);
  if (flg) { 
    PetscScalar *xx; 
    VecGetArray(x, &xx);
    VecGetOwnershipRange(x,&Istart,&Iend); 
    PetscPrintf(PETSC_COMM_WORLD, "Solution Grid (without boundary conditions):\n" );
    for (I=Istart; I<Iend; I++) {
      i = I/n;  j = I - i*n;
      PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%8.6f ", xx[I-Istart]); 
      if (j == (n-1) ) PetscSynchronizedPrintf(PETSC_COMM_WORLD,"\n");
    }
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
    VecRestoreArray(x, &xx);
  }

  PetscOptionsHasName(PETSC_NULL,"-view_sol_x",&flg);
  if (flg) { /* view solution grid in an X window */
    PetscScalar *xx;
    DA  da;
    AO  ao;
    Vec x_da; 
    DACreate2d(PETSC_COMM_WORLD,DA_NONPERIODIC,DA_STENCIL_STAR,
               n,m,PETSC_DECIDE,PETSC_DECIDE,1,0,PETSC_NULL,PETSC_NULL,&da);
    DACreateGlobalVector(da,&x_da);
    DAGetAO(da,&ao);
    VecGetOwnershipRange(x, &Istart, &Iend);
    VecGetArray(x,&xx);
    for (I=Istart; I<Iend; I++) {
      i = I; AOApplicationToPetsc(ao,1,&i);
      VecSetValues(x_da, 1, &i, &xx[I-Istart], INSERT_VALUES);
    }
    VecRestoreArray(x, &xx);
    VecAssemblyBegin(x_da); VecAssemblyEnd(x_da);
    PetscOptionsHasName(PETSC_NULL,"-view_sol_x_da",&flg);
    if (flg)  VecView(x_da,PETSC_VIEWER_STDOUT_WORLD); 
    VecView(x_da,PETSC_VIEWER_DRAW_(PETSC_COMM_WORLD)); 
    DADestroy(da); VecDestroy(x_da); 
  }

  PetscOptionsHasName(PETSC_NULL,"-view_mat_x",&flg);
  if (flg) {
    MatView(A,PETSC_VIEWER_DRAW_(PETSC_COMM_WORLD));
    /* use "-draw_pause -1" to wait after drawing until user interaction on the viewing
       window, and "-draw_pause n" to wait n seconds */ 
  }

  /* 
     Check the error
  */
  VecAXPY(x,-1,u); /* 2.3.0 */    /* VecAXPY(&neg_one,u,x) in 2.2.1 and older */
  VecNorm(x,NORM_2,&norm);
  /* Optional: Scale the norm:    norm *= sqrt(1.0/((m+1)*(n+1))); */  

  /*
     Print convergence information.  PetscPrintf() produces a single 
     print statement from all processes that share a communicator.
     An alternative is PetscFPrintf(), which prints to a file.
  */
  KSPGetIterationNumber(ksp,&its);
  PetscPrintf(PETSC_COMM_WORLD,"Norm of error %A iterations %d\n",norm,its);

  /*
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
  */
  KSPDestroy(ksp);
  VecDestroy(u);  VecDestroy(x);
  VecDestroy(b);  MatDestroy(A);

  /*
     Always call PetscFinalize() before exiting a program.  This routine
       - finalizes the PETSc libraries as well as MPI
       - provides summary and diagnostic information if certain runtime
         options are chosen (e.g., -log_summary). 
  */
  PetscFinalize();
  return 0;
}

