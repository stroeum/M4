static char help[] = "Time-dependent system of DAEs in 3d. Modified from ex13.c for illustrating how to solve DAEs. \n";

#include "MyCtx.h"
#include "MyMonitor.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
	PetscErrorCode ierr;
	MonitorCtx     usermonitor;          /* user-defined monitor context */
	AppCtx         user;                 /* user-defined work context */
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Initialize program
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	PetscInitialize(&argc,&argv,PETSC_NULL,help);
	ierr = InitCtx(&user,&usermonitor);
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Create distributed array (DMDA) to manage parallel grid and vectors
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_GHOSTED,DM_BOUNDARY_GHOSTED,DM_BOUNDARY_GHOSTED,DMDA_STENCIL_STAR,user.mx,user.my,user.mz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,19,1,PETSC_NULL,PETSC_NULL,PETSC_NULL,&user.da);CHKERRQ(ierr);
	ierr =DMSetUp(user.da); CHKERRQ(ierr);
	ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_GHOSTED,DM_BOUNDARY_GHOSTED,DM_BOUNDARY_GHOSTED,DMDA_STENCIL_STAR,user.mx,user.my,user.mz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,9,1,PETSC_NULL,PETSC_NULL,PETSC_NULL,&user.db);CHKERRQ(ierr);
	ierr =DMSetUp(user.db); CHKERRQ(ierr);

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Convert bin files into txt
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = OutputData(&user);
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Free work space.
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = DMDestroy(&user.da);CHKERRQ(ierr);
	ierr = DMDestroy(&user.db);CHKERRQ(ierr);
	free(user.x);
	free(user.y);
	free(user.z);
	
	ierr = PetscFinalize();
	PetscFunctionReturn(0);
}

