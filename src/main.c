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
	MonitorCtx     usermonitor;
	AppCtx         user;

	// For JAR RK method personalization
	PetscFunctionList TSSSPList;

	PetscInitialize(&argc,&argv,PETSC_NULL,help);
	
	char cwd[PATH_MAX];
	if (getcwd(cwd, sizeof(cwd)) != NULL) {
		PetscPrintf(PETSC_COMM_WORLD,"Current working dir: %s\n", cwd);
	} else {
		perror("getcwd() error");
		return 1;
	}

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Initialize program
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = InitCtx(&user,&usermonitor);CHKERRQ(ierr);
	if(user.chemswitch==1)
		PetscPrintf(PETSC_COMM_WORLD,"\nChemistry\t= ON\n");
	else
		PetscPrintf(PETSC_COMM_WORLD,"\nChemistry\t= OFF\n");
	if(user.collswitch==1)
		PetscPrintf(PETSC_COMM_WORLD,"Collisions\t= ON\n");
	else
		PetscPrintf(PETSC_COMM_WORLD,"Collisions\t= OFF\n");
	if(user.gradpswitch==1)
		PetscPrintf(PETSC_COMM_WORLD,"Pressure Grad\t= ON\n");
	else
		PetscPrintf(PETSC_COMM_WORLD,"Pressure Grad\t= OFF\n");
	if(user.gravswitch==1)
		PetscPrintf(PETSC_COMM_WORLD,"Gravity\t\t= ON\n\n");
	else
		PetscPrintf(PETSC_COMM_WORLD,"Gravity\t\t= OFF\n\n");
        
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Create distributed array (DMDA) to manage parallel grid and vectors
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	PetscPrintf(PETSC_COMM_WORLD,"mx = %i my = %i mz = %i\n",user.mx,user.my,user.mz);
	ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_GHOSTED,DM_BOUNDARY_GHOSTED,DM_BOUNDARY_GHOSTED,DMDA_STENCIL_STAR,user.mx,user.my,user.mz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,19,1,PETSC_NULL,PETSC_NULL,PETSC_NULL,&user.da);CHKERRQ(ierr);
	ierr = DMSetUp(user.da); CHKERRQ(ierr);
	ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_GHOSTED,DM_BOUNDARY_GHOSTED,DM_BOUNDARY_GHOSTED,DMDA_STENCIL_STAR,user.mx,user.my,user.mz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,9,1,PETSC_NULL,PETSC_NULL,PETSC_NULL,&user.db);CHKERRQ(ierr);
	ierr = DMSetUp(user.db); CHKERRQ(ierr);
	
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
	ierr = DMCreateMatrix(user.da,&J);CHKERRQ(ierr);
	ierr = TSGetSNES(ts,&snes);CHKERRQ(ierr);
		
	char version[255];
	size_t len=255;
	ierr = PetscGetVersion(version,len);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Version = %s\n",version);CHKERRQ(ierr);

	ierr = SNESSetJacobian(snes,J,J,SNESComputeJacobianDefault,(void *)&user);CHKERRQ(ierr);
	ierr = TSSetProblemType(ts,TS_NONLINEAR);CHKERRQ(ierr);
	ierr = TSSetRHSFunction(ts,r,FormFunction,(void *)&user);CHKERRQ(ierr);
	ierr = TSSetApplicationContext(ts,(void *)&user);CHKERRQ(ierr);

	ierr = TSSetMaxSteps(ts,user.maxsteps);CHKERRQ(ierr);
	ierr = TSSetMaxTime(ts,user.tf);CHKERRQ(ierr);
	ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER);CHKERRQ(ierr);
	ierr = TSMonitorSet(ts,MyTSMonitor,&user,PETSC_NULL);CHKERRQ(ierr);         // Tracks data and writes it (MyTSMonitor()) to binary at preset timesteps
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Set initial conditions
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = FormInitialSolution(u,(void *)&user);CHKERRQ(ierr);
	ierr = TSSetTime(ts,user.ti);CHKERRQ(ierr);
	ierr = TSSetTimeStep(ts,user.dt);CHKERRQ(ierr);
	ierr = TSSetSolution(ts,u);CHKERRQ(ierr);
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Set runtime options
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = TSSetFromOptions(ts);CHKERRQ(ierr);
	ierr = TSSetPostStep(ts,CFL);CHKERRQ(ierr);     // Referencing MyCtx/CFL()
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Solve nonlinear system
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = TSSolve(ts,u);CHKERRQ(ierr);
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Free work space.
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = VecDestroy(&u);CHKERRQ(ierr);
	ierr = TSDestroy(&ts);CHKERRQ(ierr);
	ierr = DMDestroy(&user.da);CHKERRQ(ierr);
	ierr = DMDestroy(&user.db);CHKERRQ(ierr);
	free(user.x);
	free(user.y);
	free(user.z);
	
	ierr = PetscFinalize();CHKERRQ(ierr);
	PetscFunctionReturn(0);
}
