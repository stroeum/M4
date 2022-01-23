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
	PetscLogDouble t0    = user->t0,t1;
	DM             da    = user->da;
	DM             db    = user->db;
	PetscReal      L     = user->L;
	PetscReal      Bt    = user->Bt;
	PetscReal      B0i   = user->B0i;
	PetscReal      B0f   = user->B0f;
	PetscReal      *x= user->x, *y=user->y, *z=user->z;
	PetscInt       i,j,k,m,xs,ys,zs,xm,ym,zm;
	PetscReal      X,Y,Z;
	PetscInt       rank;
	PetscInt       istep=user->istep, vizdstep=user->viz_dstep;
	PetscReal      tau=user->tau;
	PetscBool      xtra_out=user->xtra_out;
	Vec            V;
	Vec            localU;
	PetscReal      ****u;
	FILE           *fd;
	PetscInt       flag;
	
	PetscFunctionBegin;
	
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	// During the transition period 0->Bt, slowly introduce Bfield
	if (ptime*tau < Bt) {
		ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
		ierr = DMGetLocalVector(da,&localU );CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin(da,U,INSERT_VALUES,localU );CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd(da,U,INSERT_VALUES,localU );CHKERRQ(ierr);
		ierr = DMDAVecGetArrayDOF(da,localU ,&u );CHKERRQ(ierr);

		user->B[0] = (B0i + (B0f-B0i)*erf(ptime*tau*4/Bt - 2) + B0f)/2; // Arbitrary function which leads to a smooth transition from B0i to B0f over time Bt

		for (k=zs; k<zs+zm; k++) {
			Z = user->outZmin + z[k]*L;
			for (j=ys; j<ys+ym; j++) {
				Y = user->outYmin + y[j]*L;
				for (i=xs; i<xs+xm; i++) {
					X = user->outXmin + x[i]*L;
					//for (m=0; m<3; m++) { //Why did I have a loop here? Is this just old?
						FormInitialBField(u[k][j][i],ctx,X,Y,Z);
					//}
				}
			}
		}
	
		ierr = DMDAVecRestoreArrayDOF(da,localU ,&u );CHKERRQ(ierr);
		ierr = DMLocalToGlobalBegin(da,localU,INSERT_VALUES,U);CHKERRQ(ierr);
		ierr = DMLocalToGlobalEnd(da,localU,INSERT_VALUES,U);CHKERRQ(ierr);
		ierr = DMRestoreLocalVector(da,&localU);CHKERRQ(ierr);
	}
	
	if((istep+step)%vizdstep==0) {
		ierr = PetscTime(&t1);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Elapsed_time=%2.1e(%2.1fhr)\ttimestep %D\tt %2.3e\tdt %2.3e\n",t1-t0,(t1-t0)/3600,istep+step,ptime*tau,user->dt*tau);CHKERRQ(ierr);
        
		sprintf(fName, "%s/t.out",user->dName);
		flag = access(fName,W_OK);
		if (flag==0) fd=fopen(fName,"a");
		else         fd=fopen(fName,"w");

		if(user->isInputFile==1 && ptime==user->ti) {
			// do not record the time
		} else { 
			ierr = PetscFPrintf(PETSC_COMM_WORLD,fd,"%12.6e\n",ptime*tau); CHKERRQ(ierr);
		}
		fclose(fd);
		
		sprintf(fName, "%s/X%d.bin",user->dName,istep+step);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,fName,FILE_MODE_WRITE, &fViewer); CHKERRQ(ierr);
		ierr = VecView(U,fViewer); CHKERRQ(ierr); // try giving an integer and a long number and array figure out dim u
		ierr = PetscViewerDestroy(&fViewer); CHKERRQ(ierr);
		
		if(xtra_out) {
			ierr = DMGetGlobalVector(db,&V);CHKERRQ(ierr);
			ierr = DMGetLocalVector(da,&localU);CHKERRQ(ierr);
			
			ierr = DMGlobalToLocalBegin(da,U,INSERT_VALUES,localU);CHKERRQ(ierr);
			ierr = DMGlobalToLocalEnd(da,U,INSERT_VALUES,localU);CHKERRQ(ierr);
			ierr = DMDAVecGetArrayDOF(da,localU,&u);CHKERRQ(ierr);
			
			ierr = FormBCu(u,user);CHKERRQ(ierr);
			
			ierr = FormIntermediateFunction(u,V,user);CHKERRQ(ierr);
			
			sprintf(fName, "%s/Y%d.bin",user->dName,istep+step);
			ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,fName,FILE_MODE_WRITE, &fViewer); CHKERRQ(ierr);
			ierr = VecView(V,fViewer); CHKERRQ(ierr);
			
			ierr = DMDAVecRestoreArrayDOF(da,localU,&u);CHKERRQ(ierr);
			ierr = DMRestoreLocalVector(da,&localU);CHKERRQ(ierr);
			ierr = DMRestoreGlobalVector(db,&V);CHKERRQ(ierr);
			ierr = PetscViewerDestroy(&fViewer); CHKERRQ(ierr);
		}
	} 
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

