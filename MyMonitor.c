/*
 * MyMonitor.c
 * Created by Jeremy A Riousset on 02/04/11
 * Output solutions to binary files for with PDE solvers
 */

#include "MyMonitor.h"

#undef __FUNCT__  
#define __FUNCT__ "MyTSMonitor"
PetscErrorCode MyTSMonitor(TS ts,PetscInt step,PetscReal ptime,Vec U,void *ctx)
{
  PetscErrorCode ierr;
  char           fName[50];
  PetscViewer    fViewer;
  AppCtx         *user = (AppCtx*)ctx;
  PetscLogDouble t0 = user->t0 , t1;
  DM             db = user->db;
  PetscInt       rank;
  PetscInt       istep = user->istep, vizdstep = user->viz_dstep;
  PetscReal      tau = user->tau;
  PetscReal      dt;/*,ne;*/
  PetscBool      xtra_out = user->xtra_out;
  Vec            V;
  FILE           *fd;
  PetscInt       flag;

  PetscFunctionBegin;

  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  sprintf(fName, "%s/t.dat",user->dName);
  flag = access(fName,W_OK);
  if (flag==0) fd=fopen(fName,"a");
  else         fd=fopen(fName,"w");

  if((istep+step)%vizdstep==0) { // || (step>=9400) ) {
    ierr = TSGetTimeStep(ts,&dt);CHKERRQ(ierr);
    ierr = PetscGetTime(&t1);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Elapsed_time=%2.1e\ttimestep %D\tt %2.3e\tdt %2.3e\n",t1-t0,istep+step,ptime*tau,user->dt*tau);CHKERRQ(ierr);
    ierr = PetscFPrintf(PETSC_COMM_WORLD,fd,"%2.3e\n",ptime*tau); CHKERRQ(ierr);

    sprintf(fName, "%s/X%d.bin",user->dName,istep+step);

    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,fName,FILE_MODE_WRITE, &fViewer); CHKERRQ(ierr);
    ierr = VecView(U,fViewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fViewer); CHKERRQ(ierr);
    
    if(xtra_out) {
      ierr = DMCreateGlobalVector(db,&V);CHKERRQ(ierr);
      //PetscPrintf(PETSC_COMM_WORLD,"Extra Diagnostic Record: ...\n"); 
      ierr = FormIntermediateFunction(U,V,user);CHKERRQ(ierr);
      //PetscPrintf(PETSC_COMM_WORLD,"Extra Diagnostic Record: DONE\n"); 

      sprintf(fName, "%s/Y%d.bin",user->dName,istep+step);
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,fName,FILE_MODE_WRITE, &fViewer); CHKERRQ(ierr);
      ierr = VecView(V,fViewer); CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&fViewer); CHKERRQ(ierr);
    }
  } 
  fclose(fd);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MySNESMonitor"
/*
 MySNESMonitor - illustrate how to set user-defined monitoring routine for SNES.
 Input Parameters:
 snes - the SNES context
 its - iteration number
 fnorm - 2-norm function value (may be estimated)
 ctx - optional user-defined context for private data for the 
 monitor routine, as set by SNESMonitorSet()
 */
PetscErrorCode MySNESMonitor(SNES snes,PetscInt its,PetscReal fnorm,void *ctx)
{
  PetscErrorCode ierr;
   
  PetscFunctionBegin;
  ierr = SNESMonitorDefaultShort(snes,its,fnorm,ctx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

