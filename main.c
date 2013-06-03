/*  usage:  mpiexec -n <procs> ./main [-help] [all PETSc options] */

/* useless comment! */

static char help[] = "Nonlinear, time-dependent PDE in 3d.\n";

#include "MyCtx.h"
#include "MyMonitor.h"
#include "MyPDEs.h"
#include "petscts.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  TS             ts;  // nonlinear solver
  SNES           snes;
  Vec            u,r; // solution, residual vectors
  Mat            J;   // Jacobian matrix
  PetscErrorCode ierr;
  //MatFDColoring  matfdcoloring = PETSC_NULL;
  MonitorCtx     usermonitor;
  AppCtx         user;
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Initialize program
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscInitialize(&argc,&argv,PETSC_NULL,help);
  ierr = InitCtx(&user,&usermonitor);
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create distributed array (DMDA) to manage parallel grid and vectors
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscPrintf(PETSC_COMM_WORLD,"mx = %i my = %i mz = %i\n",user.mx,user.my,user.mz);
  ierr = DMDACreate3d(PETSC_COMM_WORLD,DMDA_BOUNDARY_GHOSTED,DMDA_BOUNDARY_GHOSTED,DMDA_BOUNDARY_GHOSTED,DMDA_STENCIL_STAR,user.mx,user.my,user.mz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,19,1,PETSC_NULL,PETSC_NULL,PETSC_NULL,&user.da);CHKERRQ(ierr);
  ierr = DMDACreate3d(PETSC_COMM_WORLD,DMDA_BOUNDARY_GHOSTED,DMDA_BOUNDARY_GHOSTED,DMDA_BOUNDARY_GHOSTED,DMDA_STENCIL_STAR,user.mx,user.my,user.mz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,9,1,PETSC_NULL,PETSC_NULL,PETSC_NULL,&user.db);CHKERRQ(ierr);

  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Extract global vectors from DMDA; then duplicate for remaining
     vectors that are the same types
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMCreateGlobalVector(user.da,&u);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&r);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create timestepping solver context
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = TSCreate(PETSC_COMM_WORLD,&ts);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSSSP);CHKERRQ(ierr);CHKERRQ(ierr);
  ierr = PetscFListAdd(&TSSSPList,"rk2jar", "TSSSPStep_RK_2_JAR", (void(*)(void))TSSSPStep_RK_2_JAR); CHKERRQ(ierr);
  ierr = PetscFListAdd(&TSSSPList,"lw", "TSSSPStep_LW", (void(*)(void))TSSSPStep_LW); CHKERRQ(ierr);
  ierr = PetscFListAdd(&TSSSPList,"lax", "TSSSPStep_LAX", (void(*)(void))TSSSPStep_LAX); CHKERRQ(ierr);
  ierr = TSSetProblemType(ts,TS_NONLINEAR);CHKERRQ(ierr);
  ierr = TSSetRHSFunction(ts,r,FormFunction,(void *)&user);CHKERRQ(ierr);
  ierr = TSSetApplicationContext(ts,(void *)&user);
#if PETSC_VERSION_LT(3,3,1)
  ierr = DMCreateMatrix(user.da,MATAIJ,&J);CHKERRQ(ierr);
#else
  ierr = DMGetMatrix(user.da,MATAIJ,&J);CHKERRQ(ierr);
#endif
  ierr = TSGetSNES(ts,&snes);CHKERRQ(ierr);
  ierr = SNESSetJacobian(snes,J,J,SNESDefaultComputeJacobian,(void *)&user);CHKERRQ(ierr);
//ierr = TSSetRHSJacobian(ts,J,J,TSDefaultComputeJacobian,(void *)&user);CHKERRQ(ierr);

  ierr = TSSetDuration(ts,user.maxsteps,user.tf);CHKERRQ(ierr);
  ierr = TSMonitorSet(ts,MyTSMonitor,&user,PETSC_NULL);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set initial conditions
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = FormInitialSolution(u,(void *)&user);CHKERRQ(ierr);
  ierr = TSSetInitialTimeStep(ts,user.ti,user.dt);CHKERRQ(ierr);
  ierr = TSSetSolution(ts,u);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set runtime options
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);
  ierr = TSSetPostStep(ts,CFL);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve nonlinear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = TSSolve(ts,u,&(user.tf));CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space.  
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
  ierr = MatDestroy(&J);CHKERRQ(ierr);
  if (matfdcoloring){
    ierr = MatFDColoringDestroy(&matfdcoloring);CHKERRQ(ierr);
  }*/
  ierr = VecDestroy(&u);CHKERRQ(ierr);     
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = DMDestroy(&user.da);CHKERRQ(ierr);
  ierr = DMDestroy(&user.db);CHKERRQ(ierr);
  free(user.x);
  free(user.y);
  free(user.z);

  ierr = PetscFinalize();
  PetscFunctionReturn(0);
}
