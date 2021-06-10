/*
 * MyCtx.c
 * Created by Jeremy A Riousset on 02/04/11
 * User defined contexts initialziation for use with DAE solvers
 */

#include "MyCtx.h"

#undef __FUNCT__
#define __FUNCT__ "TSSSPStep_RK_2_JAR"
/*
 * TSSSPRKS2 - Optimal second order SSP Runge-Kutta method, low-storage, c_eff=(s-1)/s
 */
PetscErrorCode TSSSPStep_RK_2_JAR(TS ts,PetscReal t0,PetscReal dt,Vec sol)
{
	AppCtx         *user;
	PetscErrorCode ierr;
	Vec            W,F;
	PetscBool      smoothing;
	
	PetscFunctionBegin;
	
	ierr = TSGetApplicationContext(ts,(void *)&user);CHKERRQ(ierr);
	smoothing = user->smoothing;
	
	ierr = VecDuplicate(sol,&W);CHKERRQ(ierr);
	ierr = VecDuplicate(sol,&F);CHKERRQ(ierr);
	ierr = VecCopy(sol,W);CHKERRQ(ierr);
	
	ierr = TSComputeRHSFunction(ts,t0,W,F);CHKERRQ(ierr);
	ierr = VecAXPY(W,dt/2.0,F);CHKERRQ(ierr);
	
	ierr = TSComputeRHSFunction(ts,t0+dt/2.0,W,F);CHKERRQ(ierr);
	if (smoothing) {ierr = WeightedAverage(F,user);CHKERRQ(ierr);}
	ierr = VecAXPY(sol,dt,F);CHKERRQ(ierr);
	ierr = VecDestroy(&W); CHKERRQ(ierr);
	ierr = VecDestroy(&F); CHKERRQ(ierr);
	user->cnt++;
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TSSSPStep_LAX"
/*MC
 TSSSPLAX - First order Lax method
 M*/
PetscErrorCode TSSSPStep_LAX(TS ts,PetscReal t0,PetscReal dt,Vec sol)
{
	AppCtx         *user;
	PetscErrorCode ierr;
	Vec            W,F;
	
	PetscFunctionBegin;
	
	ierr = TSGetApplicationContext(ts,(void *)&user);CHKERRQ(ierr);
	
	ierr = VecDuplicate(sol,&W);CHKERRQ(ierr);
	ierr = VecDuplicate(sol,&F);CHKERRQ(ierr);
	ierr = VecCopy(sol,W);CHKERRQ(ierr);
	ierr = TSComputeRHSFunction(ts,t0,W,F);CHKERRQ(ierr);
	
	ierr = Average(W,user); CHKERRQ(ierr);
	ierr = VecWAXPY(sol,dt,F,W);CHKERRQ(ierr);
	
	ierr = VecDestroy(&W); CHKERRQ(ierr);
	ierr = VecDestroy(&F); CHKERRQ(ierr);
	user->cnt++;
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TSSSPStep_LW"
/*MC
 TSSSPLW - Second order Lax-Wendroff method
 M*/
PetscErrorCode TSSSPStep_LW(TS ts,PetscReal t0,PetscReal dt,Vec sol)
{
	AppCtx         *user;
	PetscErrorCode ierr;
	Vec            W,F;
	
	PetscFunctionBegin;
	
	ierr = TSGetApplicationContext(ts,(void *)&user);CHKERRQ(ierr);
	
	ierr = VecDuplicate(sol,&W);CHKERRQ(ierr);
	ierr = VecDuplicate(sol,&F);CHKERRQ(ierr);
	ierr = VecCopy(sol,W);CHKERRQ(ierr);
	
	ierr = TSComputeRHSFunction(ts,t0,W,F);CHKERRQ(ierr); // W^n = U^n -> F = f(U^n)
	ierr = Average(W,user); CHKERRQ(ierr); // W^n = <U^n>
	ierr = VecAXPY(W,dt/2.0,F);CHKERRQ(ierr); // W^n = <U^n> + dt/2*f(U^n) = U^{n+1/2}
	
	ierr = TSComputeRHSFunction(ts,t0+dt/2.0,W,F);CHKERRQ(ierr); // F = f(U^{n+1/2})
	ierr = VecAXPY(sol,dt,F);CHKERRQ(ierr); // sol = sol + dt*F = U^n + dt f(U^{n+1/2})
	ierr = VecDestroy(&W); CHKERRQ(ierr);
	ierr = VecDestroy(&F); CHKERRQ(ierr);
	user->cnt++;
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "InitCtx"
/* ------------------------------------------------------------------- */
/* 
 InitCtx - Assign values to user-defined context.
 
 Input Parameters:
 .  user - user-defined context, as set by SNESSetFunction()
 .  usermonitor - user-defined monitor context
 
 Output Parameter:
 */
PetscErrorCode InitCtx(AppCtx *user, MonitorCtx *usrmnt)
{
	PetscErrorCode ierr;
	PetscInt       i=0, j=0, k=0, rank=0;
	PetscInt       mx_lo,my_lo,mz_lo;
	PetscInt       mx_in,my_in,mz_in;
	PetscInt       mx_hi,my_hi,mz_hi;
	PetscReal      pct=1.1;
	PetscReal      dx=0, dy=0, dz=0;
	PetscReal      X=0, Y=0, Z=0;
	svi            s = {// Index for variables with no d*/dt (svi)
		{0,1,2}, // index of vex,vey,vez
		{3,4,5}, // index of Jx,Jy,Jz
		{6,7,8}  // index of Ex,Ey,Ez
	};
	dvi            d = {// Index for variables with d*/dt (dvi)
		{0,1,2},    // index of ni[O2p],ni[CO2p],ni[Op]
		{
		{3,4,5},  // index of x-,y-,z-components of vi[O2p]
		{6,7,8},  // index of x-,y-,z-components of vi[CO2p]
		{9,10,11},// index of x-,y-,z-components of vi[Op]
		},
		{12,13,14}, // index of pi[O2p],pi[CO2p],pi[Op]
		15,         // index of pe
		{16,17,18}  // index of Bx,By,Bz
	};
	
	char           InputFile[PETSC_MAX_PATH_LEN];
	char           *pdot;
	char           dot = '.';
	char           *pstep = NULL;
	
	user->s = s;
	user->d = d;
	user->da = PETSC_NULL;
	user->db = PETSC_NULL;
	user->cnt = 0;
	
	user->gama[O2p]  = 1.0 + 2.0/6.0;
	user->gama[CO2p] = 1.0 + 2.0/9.0;
	user->gama[Op]   = 1.0 + 2.0/3.0;
	user->gama[e]    = 1.0 + 2.0/3.0;
	
	PetscFunctionBegin;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	
	ierr = PetscTime(&user->t0);CHKERRQ(ierr);
	user->ldName = sprintf(user->dName,"output/");
	ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-output_folder",user->dName,PETSC_MAX_PATH_LEN-1,PETSC_NULL);CHKERRQ(ierr);
	
	user->isInputFile = PETSC_FALSE;
	user->ldName = sprintf(user->InputFile,"X0.bin");
	ierr = PetscOptionsGetString(PETSC_NULL,PETSC_NULL,"-input_file",user->InputFile,PETSC_MAX_PATH_LEN-1,&user->isInputFile);CHKERRQ(ierr);
	
	usrmnt->drawcontours = PETSC_FALSE;
	ierr = PetscOptionsHasName(PETSC_NULL,PETSC_NULL,"-drawcontours",&usrmnt->drawcontours);CHKERRQ(ierr);
	user->maxsteps  = 100;
	ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-ts_max_steps",&user->maxsteps,PETSC_NULL);CHKERRQ(ierr);
	user->bcType = 0;
	ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-bc_type",&user->bcType,PETSC_NULL);CHKERRQ(ierr);
	user->bcRec = PETSC_TRUE;
	
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-mi0",&user->mi[0],PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-mi1",&user->mi[1],PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-mi2",&user->mi[2],PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-me",&user->me,PETSC_NULL);CHKERRQ(ierr);
	
	user->un[0] = 0.0;
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-WindNS",&user->un[0],PETSC_NULL);CHKERRQ(ierr);
	user->un[1] = 0.0;
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-WindEW",&user->un[1],PETSC_NULL);CHKERRQ(ierr);
	user->un[2] = 0.0;
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-WindUD",&user->un[2],PETSC_NULL);CHKERRQ(ierr);
	
	user->ui[0] = 0.0;
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-iWindNS",&user->ui[0],PETSC_NULL);CHKERRQ(ierr);
	user->ui[1] = 0.0;
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-iWindEW",&user->ui[1],PETSC_NULL);CHKERRQ(ierr);
	user->ui[2] = 0.0;
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-iWindUD",&user->ui[2],PETSC_NULL);CHKERRQ(ierr);
	
	user->xtra_out = PETSC_FALSE;
	ierr = PetscOptionsHasName(PETSC_NULL,PETSC_NULL,"-ExtraDiagnostics",&user->xtra_out);CHKERRQ(ierr);
	user->smoothing = PETSC_FALSE;
	ierr = PetscOptionsHasName(PETSC_NULL,PETSC_NULL,"-smoothing",&user->smoothing);CHKERRQ(ierr);
	user->limiters = PETSC_FALSE;
	ierr = PetscOptionsHasName(PETSC_NULL,PETSC_NULL,"-limiters",&user->limiters);CHKERRQ(ierr);
	user->blockers = PETSC_FALSE;
	ierr = PetscOptionsHasName(PETSC_NULL,PETSC_NULL,"-blockers",&user->blockers);CHKERRQ(ierr);
	user->BfieldType = 0;
	ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-B_field_type",&user->BfieldType,PETSC_NULL);CHKERRQ(ierr);
	user->B[0] = 0.0;
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-B0",&user->B[0],PETSC_NULL);CHKERRQ(ierr);
	user->B[1] = 0.0;
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-B1",&user->B[1],PETSC_NULL);CHKERRQ(ierr);
	user->B[2] = 0.0;
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-B2",&user->B[2],PETSC_NULL);CHKERRQ(ierr);
	user->B[3] = 0.0;
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-B3",&user->B[3],PETSC_NULL);CHKERRQ(ierr);
	
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-Rmars",&user->rM,PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-Mmars",&user->mM,PETSC_NULL);CHKERRQ(ierr);
	
	user->inXmin = 0.0;
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-in_x_min",&user->inXmin,PETSC_NULL);CHKERRQ(ierr);
	user->inXmax = 1.0;
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-in_x_max",&user->inXmax,PETSC_NULL);CHKERRQ(ierr);
	user->inYmin = 0.0;
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-in_y_min",&user->inYmin,PETSC_NULL);CHKERRQ(ierr);
	user->inYmax = 1.0;
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-in_y_max",&user->inYmax,PETSC_NULL);CHKERRQ(ierr);
	user->inZmin = 0.0;
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-in_z_min",&user->inZmin,PETSC_NULL);CHKERRQ(ierr);
	if(user->inZmin<0) user->inZmin = 0.0;
	user->inZmax = 1.0;
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-in_z_max",&user->inZmax,PETSC_NULL);CHKERRQ(ierr);
	user->outXmin = user->inXmin;
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-out_x_min",&user->outXmin,PETSC_NULL);CHKERRQ(ierr);
	user->outXmax = user->inXmax;
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-out_x_max",&user->outXmax,PETSC_NULL);CHKERRQ(ierr);
	user->outYmin = user->inYmin;
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-out_y_min",&user->outYmin,PETSC_NULL);CHKERRQ(ierr);
	user->outYmax = user->inYmax;
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-out_y_max",&user->outYmax,PETSC_NULL);CHKERRQ(ierr);
	user->outZmin = user->inZmin;
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-out_z_min",&user->outZmin,PETSC_NULL);CHKERRQ(ierr);
	if(user->outZmin<0) user->outZmin = 0.0;
	user->outZmax = user->inZmax;
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-out_z_max",&user->outZmax,PETSC_NULL);CHKERRQ(ierr);
	
	user->vDamping = PETSC_FALSE;
	ierr = PetscOptionsHasName(PETSC_NULL,PETSC_NULL,"-vDamping",&user->vDamping);CHKERRQ(ierr);
	user->ZL    = user->inZmin;
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-Damping_lower_alt",&user->ZL,PETSC_NULL);CHKERRQ(ierr);
	user->ZU    = user->inZmax;
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-Damping_upper_alt",&user->ZU,PETSC_NULL);CHKERRQ(ierr);
	user->lambda= 1.0;
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-Damping_scaleheight",&user->lambda,PETSC_NULL);CHKERRQ(ierr);
	
	if(user->outXmin >= user->inXmin) user->outXmin = user->inXmin;
	if(user->outXmax <= user->inXmax) user->outXmax = user->inXmax;
	if(user->outYmin >= user->inYmin) user->outYmin = user->inYmin;
	if(user->outYmax <= user->inYmax) user->outYmax = user->inYmax;
	if(user->outZmin >= user->inZmin) user->outZmin = user->inZmin;
	if(user->outZmax <= user->inZmax) user->outZmax = user->inZmax;
	
	/* Commented by Kellen to remove vizbox
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-vizbox_x_min",&user->vizbox[0],PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-vizbox_x_max",&user->vizbox[1],PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-vizbox_y_min",&user->vizbox[2],PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-vizbox_y_max",&user->vizbox[3],PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-vizbox_z_min",&user->vizbox[4],PETSC_NULL);CHKERRQ(ierr);
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-vizbox_z_max",&user->vizbox[5],PETSC_NULL);CHKERRQ(ierr);
	*/
	
	user->viz_dstep = 0;
	ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-viz_dstep",&user->viz_dstep,PETSC_NULL);CHKERRQ(ierr);
	
	mx_in = 2;
	ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-da_grid_x",&mx_in,PETSC_NULL);CHKERRQ(ierr);
	my_in = 2;
	ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-da_grid_y",&my_in,PETSC_NULL);CHKERRQ(ierr);
	mz_in = 2; 
	ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-da_grid_z",&mz_in,PETSC_NULL);CHKERRQ(ierr);
	
	user->ti = 0.0; 
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-ts_init_time",&user->ti,PETSC_NULL);CHKERRQ(ierr);
	user->tf = 1.0;
	ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-ts_max_time",&user->tf,PETSC_NULL);CHKERRQ(ierr);

    user->chemswitch = 0;
    ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-chem_switch",&user->chemswitch,PETSC_NULL);CHKERRQ(ierr);
    user->collswitch = 0;
    ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-coll_switch",&user->collswitch,PETSC_NULL);CHKERRQ(ierr);
    user->gravswitch = 0;
    ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-grav_switch",&user->gravswitch,PETSC_NULL);CHKERRQ(ierr);
    
    user->dt = 1.0e-9;
    ierr = PetscOptionsGetReal(PETSC_NULL,PETSC_NULL,"-dt",&user->dt,PETSC_NULL);CHKERRQ(ierr);
        
	if(user->viz_dstep<=0 || (PetscInt)(user->viz_dstep)!=user->viz_dstep){
		PetscPrintf(PETSC_COMM_WORLD,"ERROR: Invalid value of viz_dstep: viz_dstep=1");
		user->viz_dstep = 1;
	}
	
	user->gM = G*user->mM/(user->rM*user->rM);
	
	user->Lx = user->inXmax-user->inXmin;
	user->Ly = user->inYmax-user->inYmin;
	user->Lz = user->inZmax-user->inZmin;
	user->L  = 1.0e3;
	
	dx = user->Lx/(PetscReal)(mx_in-1);
	dy = user->Ly/(PetscReal)(my_in-1);
	dz = user->Lz/(PetscReal)(mz_in-1);
	
	mx_lo = (PetscInt)ceil(log((pct-1.0)*(user->inXmin-user->outXmin)/(pct*dx) + 1)/log(pct)); 
	my_lo = (PetscInt)ceil(log((pct-1.0)*(user->inYmin-user->outYmin)/(pct*dy) + 1)/log(pct));
	mz_lo = (PetscInt)ceil(log((pct-1.0)*(user->inZmin-user->outZmin)/(pct*dz) + 1)/log(pct));
	
	mx_hi = (PetscInt)floor(log((pct-1.0)*(user->outXmax-user->inXmax)/(pct*dx) + 1)/log(pct)); 
	my_hi = (PetscInt)floor(log((pct-1.0)*(user->outYmax-user->inYmax)/(pct*dy) + 1)/log(pct));
	mz_hi = (PetscInt)floor(log((pct-1.0)*(user->outZmax-user->inZmax)/(pct*dz) + 1)/log(pct));
	
	// Calculate dimension of the total system from outerbound to outerbound
	X = user->inXmin - pct*dx * (pow(pct,mx_lo)-1) / (pct-1);
	if (X > user->outXmin) mx_lo++;
	X = user->inXmax + pct*dx * (pow(pct,mx_hi)-1) / (pct-1);
	if (X < user->outXmax) mx_hi++;
	Y = user->inYmin - pct*dy * (pow(pct,my_lo)-1) / (pct-1);
	if (Y > user->outYmin) my_lo++;
	Y = user->inYmax + pct*dy * (pow(pct,my_hi)-1) / (pct-1);
	if (Y < user->outYmax) my_hi++;
	Z = user->inZmin - pct*dz * (pow(pct,mz_lo)-1) / (pct-1);
	if (Z > user->outZmin) mz_lo++;
	Z = user->inZmin - pct*dz * (pow(pct,mz_lo)-1) / (pct-1);
	if (Z < 0.0) mz_lo--;
	Z = user->inZmax + pct*dz * (pow(pct,mz_hi)-1) / (pct-1);
	if (Z < user->outZmax) mz_hi++;
	
	user->mx = mx_lo + mx_in + mx_hi;
	user->my = my_lo + my_in + my_hi;
	user->mz = mz_lo + mz_in + mz_hi;
	
	user->x = malloc(user->mx*sizeof(PetscReal));
	user->y = malloc(user->my*sizeof(PetscReal));
	user->z = malloc(user->mz*sizeof(PetscReal));
	
	// Calculate the position
	user->outXmin = user->inXmin - pct*dx * (pow(pct,mx_lo)-1) / (pct-1);
	for (i=0; i<user->mx; i++) {
		if (i < mx_lo) X = user->inXmin - pct*dx * (pow(pct,mx_lo-i)-1) / (pct-1);
		else if (i >= mx_lo && i < mx_lo+mx_in) X = user->inXmin + (i-mx_lo)*dx;
		else X = user->inXmax + pct*dx * (pow(pct,i-mx_lo-mx_in+1)-1) / (pct-1);
		user->x[i] = (X - user->outXmin)/user->L;
	}
	user->outXmax = X;
	
	user->outYmin = user->inYmin - pct*dy * (pow(pct,my_lo)-1) / (pct-1);
	for (j=0; j<user->my; j++) {
		if (j < my_lo) Y = user->inYmin - pct*dy * (pow(pct,my_lo-j)-1) / (pct-1);
		else if (j >= my_lo && j < my_lo+my_in) Y = user->inYmin + (j-my_lo)*dy;
		else Y = user->inYmax + pct*dy * (pow(pct,j-my_lo-my_in+1)-1) / (pct-1);
		user->y[j] = (Y - user->outYmin)/user->L;
	}
	user->outYmax = Y;
	
	user->outZmin = user->inZmin - pct*dz * (pow(pct,mz_lo)-1) / (pct-1);
	for (k=0; k<user->mz; k++) {
		if (k < mz_lo) Z = user->inZmin - pct*dz * (pow(pct,mz_lo-k)-1) / (pct-1);
		else if (k >= mz_lo && k < mz_lo+mz_in) Z = user->inZmin + (k-mz_lo)*dz;
		else Z = user->inZmax + pct*dz * (pow(pct,k-mz_lo-mz_in+1)-1) / (pct-1);
		user->z[k] = (Z - user->outZmin)/user->L;
	}
	user->outZmax = Z;
	
	user->tau = PetscSqrtScalar(user->L/user->gM);
	user->n0  = user->me/(qe*qe*user->L*user->L*mu0);
	user->v0  = PetscSqrtScalar(user->gM*user->L);
	user->p0  = (user->me*user->me*user->gM)/(qe*qe*user->L*mu0);
	user->T0  = user->p0/(kB*user->n0 );
	user->B0  = user->me/(qe*user->tau);
	
	user->dt    = user->dt/user->tau;     // Initial time step (need to write this in main.in)
	user->eps.n = 1.0e3/user->n0; // min density allowed in the domain: 10^-3 cm-3 = 1e3 m^-3
	user->eps.p = 1e-15/user->p0; // min pressure allowed in the domain
	user->eps.v = .01*299742458/user->v0; // MAX velocity allowed in the domain: 1% of the speed of light
	
	// Retrieve information from initial conditions
	user->istep = 0;
	if (user->isInputFile) {
		/*
		 * Retrieve step number from the namefile of Xnnnnnnn.dat
		 * 1. Find the position of the "." character
		 * 2. Define the length of nnnnnnn
		 * 3. Create a pointer on char pstep and point it on a string from the second character up to "." character of Xnnnnnn.dat
		 * 4. Add end of string character to the pstep
		 * 5. Convert to integer and store initial step number into user->istep
		 */
		sprintf(InputFile,"%s",user->InputFile);
		pdot = strchr(InputFile,dot); // pointer on the "." character
		if(pdot) {
			size_t len = pdot - InputFile -1; // string length
			pstep = malloc(len);
			memcpy(pstep, &InputFile[1],len);
			pstep[len] = '\0';
			user->istep = atoi(pstep);
			
			/*
			 * Retrieve initial time from step number and file t.out containing all the time sequence.
			 * 1. Retrieve the line number
			 * 2. Scan the file t.out until this line
			 * 3. Store the time in the variable t
			 * 4. Store the initial time in user->ti
			 */
			PetscInt  tg_line;
			PetscInt  line=1;
			tg_line = 1 + user->istep/user->viz_dstep;
			FILE *pFile;
			float t=0;
			pFile = fopen("output/t.out","r");
			while(pFile && line<=tg_line) {
				fscanf(pFile,"%e",&t);
				line++;
			}
			line--;
			fclose(pFile);
			user->ti = t/user->tau; 
			
			PetscPrintf(PETSC_COMM_WORLD,"Resuming simulation from:\n * File: %s\n * Step: %d\n * Time: %e s\n * Line: %d in output/t.out\n", user->InputFile, user->istep, user->ti, line);

		} else {
			PetscPrintf(PETSC_COMM_WORLD,"Unreadable input file\n");
			exit(1);
		}
	}
	
	/* Commented by Kellen to remove vizbox
	if(user->vizbox[0]<user->outXmin || user->vizbox[1]>user->outXmax || user->vizbox[2]<user->outYmin || user->vizbox[3]>user->outYmax || user->vizbox[4]<user->outZmin || user->vizbox[5]>user->outZmax){
		PetscPrintf(PETSC_COMM_WORLD,"WARNING: Visualization area greater than simulation domain:\n---> the vizualization box is adjusted to fit the maximum dimension of the simulation domain.\n");
		if (user->vizbox[0]<user->outXmin) user->vizbox[0] = user->outXmin;
		if (user->vizbox[1]>user->outXmax) user->vizbox[1] = user->outXmax;
		if (user->vizbox[2]<user->outYmin) user->vizbox[2] = user->outYmin;
		if (user->vizbox[3]>user->outYmax) user->vizbox[3] = user->outYmax;
		if (user->vizbox[4]<user->outZmin) user->vizbox[4] = user->outZmin;
		if (user->vizbox[5]>user->outZmax) user->vizbox[5] = user->outZmax;
	}
	*/
	
	// Read input files
	if(rank==0) {
		//Changed by Kellen so that the program is run from parent directory M4, not M4/bin
		//ReadTable(user->RefProf,Np_REF,"../input/Profiles.in");
		//ReadTable(user->RefPart,4,"../input/Partition.in");
		ReadTable(user->RefProf,Np_REF,"./input/Profiles.in");
		ReadTable(user->RefPart,4,"./input/Partition.in");
	}
	MPI_Bcast(&(user->RefProf[0][0]),Np_REF*Nz_REF, MPIU_REAL, 0, PETSC_COMM_WORLD);
	MPI_Bcast(&(user->RefPart[0][0]),     4*Nz_REF, MPIU_REAL, 0, PETSC_COMM_WORLD);
	
	// Display diagnostics
	PetscPrintf(PETSC_COMM_WORLD,"\nmx = [%i %i %i]:%i\nmy = [%i %i %i]:%i\nmz = [%i %i %i]:%i\n\n",mx_lo,mx_in,mx_hi,user->mx, my_lo,my_in,my_hi,user->my, mz_lo,mz_in,mz_hi,user->mz); 
	PetscPrintf(PETSC_COMM_WORLD,"      O2+\t\tCO2+\t\tO+\t\te\n");
	PetscPrintf(PETSC_COMM_WORLD,"m   = [%12.6e\t%12.6e\t%12.6e\t%12.6e]\n",user->mi[0],user->mi[1],user->mi[2],user->me);
	PetscPrintf(PETSC_COMM_WORLD,"v0  = [%12.6e]\n",user->v0);
	PetscPrintf(PETSC_COMM_WORLD,"tau = [%12.6e]\n",user->tau);
	PetscPrintf(PETSC_COMM_WORLD,"n0  = [%12.6e]\n",user->n0);
	PetscPrintf(PETSC_COMM_WORLD,"p0  = [%12.6e]\n",user->p0);
	PetscPrintf(PETSC_COMM_WORLD,"B0  = [%12.6e]\n",user->B0);
	PetscPrintf(PETSC_COMM_WORLD,"E0  = [%12.6e]\n",user->v0*user->B0);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "OutputData"
PetscErrorCode OutputData(void* ptr)
{
	PetscErrorCode ierr;
	AppCtx         *user    = (AppCtx*)ptr;
	PetscInt       maxsteps = user->maxsteps;
	PetscInt       vizdstep = user->viz_dstep;
	PetscBool      xtra_out = user->xtra_out;
	PetscReal      Xmin     = user->outXmin, Ymin = user->outYmin, Zmin = user->outZmin;
	PetscReal      n0       = user->n0, v0 = user->v0, p0 = user->p0, B0 = user->B0;
	PetscReal      *x       = user->x, *y = user->y, *z = user->z;
	PetscReal      DX,DY,DZ,DV;
	PetscInt       mx       = user->mx, my = user->my, mz = user->mz;
	DM             da       = user->da;
	DM             db       = user->db;
	PetscReal      L        = user->L;
	float          t;
	dvi            d        = user->d;
	svi            s        = user->s;
	/* Commented by Kellen to remove vizbox
	PetscReal      vizbox[6]= {
		user->vizbox[0],user->vizbox[1],
		user->vizbox[2],user->vizbox[3],
		user->vizbox[4],user->vizbox[5]
	};
	*/
	
	Vec            U        = NULL;
	Vec            V        = NULL;
	PetscReal      ****u,****v = NULL;
	FILE           *pFile,*nFile,*tFile,*NFile;
	char           fName[PETSC_MAX_PATH_LEN];
	PetscViewer    fViewer;
	PetscInt       rank,step,flag;
	PetscInt       i,j,k,l,m,xs,ys,zs,xm,ym,zm;
	PetscInt       imin,imax,jmin,jmax,kmin,kmax; 
	PetscReal      ne = 0,ni[3],ve[3],vi[3][3],pe,pi[3],J[3],E[3],B[3]; 
	PetscReal      Ne,Ni[3],Fe[6],Fi[3][6];
	PetscReal      X,Y,Z;
	
	PetscFunctionBegin;
	
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	
	/* Commented by Kellen to remove vizbox
	PetscPrintf(PETSC_COMM_WORLD,"box: [%2.1f,%2.1f,%2.1f] to [%2.1f,%2.1f,%2.1f] (km)\n",vizbox[0]/1e3,vizbox[2]/1e3,vizbox[4]/1e3,vizbox[1]/1e3,vizbox[3]/1e3,vizbox[5]/1e3);
	*/
	
	// Store total number of charge carriers in a file named diagnotics.dat
	int kTop=0;
	Z = Zmin + z[kTop]*L;
	while(Z<300e3) {
		kTop++;
		Z = Zmin + z[kTop]*L;
	}
	PetscPrintf(PETSC_COMM_WORLD,"kTop=%d, fluxAlt=%f\n", kTop, Z*1e-3);
	
	// Commented by Kellen to remove .dat
	/*
	sprintf(fName,"%s/diagnostics.dat",user->vName);
	flag  = access(fName,W_OK);
	if (flag==0) nFile = fopen(fName,"a");
	else         nFile = fopen(fName,"w");
	ierr = PetscFPrintf(PETSC_COMM_WORLD,nFile,"t           \tN(O2+)      \tF.w(O2+)    \tF.e(O2+)    \tF.s(O2+)    \tF.n(O2+)    \tF.b(O2+)    \tF.t(O2+)    \tN(CO2+)     \tF.w(CO2+)   \tF.e(CO2+)   \tF.s(CO2+)   \tF.n(CO2+)   \tF.b(CO2+)   \tF.t(CO2+)   \tN(O+)       \tF.w(O+)     \tF.e(O+)     \tF.s(O+)     \tF.n(O+)     \tF.b(O+)     \tF.t(O+)     \tN(e)        \tF.w(e)      \tF.e(e)      \tF.s(e)      \tF.n(e)      \tF.b(e)      \tF.t(e)      \n");CHKERRQ(ierr);
	
	sprintf(fName,"%s/densities.dat",user->vName);
	flag  = access(fName,W_OK);
	if (flag==0) NFile = fopen(fName,"a");
	else         NFile = fopen(fName,"w");
	ierr = PetscFPrintf(PETSC_COMM_WORLD,NFile,"Z (m)       \tn.O2+  (m-3)\tn.CO2+ (m-3)\tn.O+   (m-3)\tn.e    (m-3)\n");CHKERRQ(ierr);
	*/
	
	// Changed by Kellen to run from root directory not bin
	//tFile = fopen("../output/t.out","r");
	tFile = fopen("./output/t.out","r");
	
	// Get local grid boundaries
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	for (step=0; step<=maxsteps; step++) {
		
		if(step%vizdstep==0) {
			
			ierr = PetscPrintf(PETSC_COMM_WORLD,"timestep %D\n",step);CHKERRQ(ierr);
			
			// Retrieve current time
			fscanf(tFile,"%e",&t);
			
			// Read binary file containing the solution
			sprintf(fName,"%s/X%d.bin",user->dName,step);
			flag  = access(fName,F_OK);
			if (flag != 0) break;
			ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,fName,FILE_MODE_READ, &fViewer); CHKERRQ(ierr);
			ierr = VecCreate(PETSC_COMM_WORLD,&U);CHKERRQ(ierr);
			ierr = VecSetBlockSize(U,19);CHKERRQ(ierr);
			ierr = VecLoad(U,fViewer); CHKERRQ(ierr);
			ierr = PetscViewerDestroy(&fViewer);CHKERRQ(ierr);
			if(xtra_out) {
				sprintf(fName, "%s/Y%d.bin",user->dName,step);
				flag  = access(fName,F_OK);
				if (flag != 0) break;
				ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,fName,FILE_MODE_READ, &fViewer); CHKERRQ(ierr);
				ierr = VecCreate(PETSC_COMM_WORLD,&V); CHKERRQ(ierr);
				ierr = VecSetBlockSize(V,9);CHKERRQ(ierr);
				ierr = VecLoad(V,fViewer); CHKERRQ(ierr);
				ierr = PetscViewerDestroy(&fViewer);CHKERRQ(ierr);
			}
			
			// Get pointers to vector data
			//ierr = DMDAVecGetArrayDOF(da,U,&u);CHKERRQ(ierr);
			//if(xtra_out) {ierr = DMDAVecGetArrayDOF(db,V,&v);CHKERRQ(ierr);}

			//ierr = DMDAVecRestoreArrayDOF(da,U,&u);CHKERRQ(ierr);
			//ierr = VecDestroy(&U);CHKERRQ(ierr);
			//if(xtra_out) {
				//ierr = DMDAVecRestoreArrayDOF(db,V,&v);CHKERRQ(ierr);
				//ierr = VecDestroy(&V);CHKERRQ(ierr);
			//}
		}
	}
	fclose(NFile);
	fclose(tFile);
	fclose(nFile);
	PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
/* Calculate the Courant-Friedrich-Lewy condition                      */
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "CFL"
PetscErrorCode CFL(TS ts)
{
	PetscErrorCode ierr;
	AppCtx         *user;
	
	PetscFunctionBegin;
	ierr = TSGetApplicationContext(ts,(void *)&user);CHKERRQ(ierr);
	ierr = TSSetTimeStep(ts,user->dt);CHKERRQ(ierr);
	if (user->dt*user->tau < 1e-12) { 
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Time step is < 1 ps.\n");CHKERRQ(ierr);
		ierr = TSSetConvergedReason(ts,TS_CONVERGED_USER);CHKERRQ(ierr);
	};
	PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
/* Calculate the magnetic field of multiple vertical dipoles           */
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "ARCADES"
/*
 * Arcades with the following configuration
 * X X X
 * 0 0 0
 */
PetscReal Arcades(PetscReal x, PetscReal y, PetscReal z, PetscInt m)
{
	PetscInt  i,j,k,s;
	PetscReal B = 0., mu = 1e16;
	PetscReal xs[3] = {-100e3,0e3,100e3};
	PetscReal ys[2] = {-50e3,50e3};
	PetscReal zs[1] = {-20e3};
	PetscInt Ms=sizeof(xs)/sizeof(xs[0]), Ns=sizeof(ys)/sizeof(ys[0]), Ps=sizeof(zs)/sizeof(zs[0]);
	B = 0;
	for (i=0; i<Ms; i++) {
		for (j=0; j<Ns; j++) {
			for (k=0; k<Ps; k++) {
				s = j; 
				B += pow(-1,s%2) * V_Dipole(mu,xs[i],ys[j],zs[k],x,y,z,m);
			}
		}
	}
	return B;
}

/* ------------------------------------------------------------------- */
/* Calculate the magnetic field of multiple Arcades                    */
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "MULTIARCADES"
/*
 * MultipleArcades
 */
PetscReal MultiArcades(PetscReal x, PetscReal y, PetscReal z, PetscInt m)
{
	PetscInt  i,j,k,s;
	PetscReal B = 0., mu = 1e16;
	PetscReal xs[3] = {-100e3,0e3,100e3};
	PetscReal ys[3] = {-250e3,0e3,250e3};
	PetscReal zs[1] = {-20e3};
	PetscInt Ms=sizeof(xs)/sizeof(xs[0]), Ns=sizeof(ys)/sizeof(ys[0]), Ps=sizeof(zs)/sizeof(zs[0]);
	B = 0;
	for (i=0; i<Ms; i++) {
		for (j=0; j<Ns; j++) {
			for (k=0; k<Ps; k++) {
				s = j; 
				B += pow(-1,s%2) * V_Dipole(mu,xs[i],ys[j],zs[k],x,y,z,m);
			}
		}
	}
	return B;
}

/* ------------------------------------------------------------------- */
/* Calculate the magnetic field of a vertical dipole                   */
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "V_DIPOLE"
/*
 * Vertical dipole
 */
PetscReal V_Dipole(PetscReal mu, PetscReal xs, PetscReal ys, PetscReal zs, PetscReal x, PetscReal y, PetscReal z, PetscInt m)
{
	PetscReal R;
	PetscReal X,Y,Z;
	X = x-xs;
	Y = y-ys;
	Z = z-zs;
	R = PetscSqrtScalar(X*X+Y*Y+Z*Z);
	
	if (m==0) return mu0/(4*M_PI) * mu * pow(1/R,3) * 3.0*X*Z/(R*R) ;
	if (m==1) return mu0/(4*M_PI) * mu * pow(1/R,3) * 3.0*Y*Z/(R*R) ;
	if (m==2) return mu0/(4*M_PI) * mu * pow(1/R,3) *(3.0*Z*Z/(R*R) - 1.0);
	return 0;
}

/* ------------------------------------------------------------------- */
/* Calculate the magnetic field of a horizontal dipole                 */
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "H_DIPOLE"
/*
 * Horizontal dipole
 */
PetscReal H_Dipole(PetscReal mu, PetscReal xs, PetscReal ys, PetscReal zs, PetscReal x, PetscReal y, PetscReal z, PetscInt m)
{
	PetscReal R;
	PetscReal X,Y,Z;
	X = x-xs;
	Y = y-ys;
	Z = z-zs;
	R = PetscSqrtScalar(X*X+Y*Y+Z*Z);
	if (m==0) return  mu0/(4*M_PI) * mu * pow(1/R,3) *(3.0*Z*Z/(R*R) - 1.0);
	if (m==1) return  mu0/(4*M_PI) * mu * pow(1/R,3) * 3.0*Y*Z/(R*R) ;
	if (m==2) return -mu0/(4*M_PI) * mu * pow(1/R,3) * 3.0*X*Z/(R*R) ;
	return 0;
}

/* ------------------------------------------------------------------- */
/* Calculate the vectorial (cross) product P =  M x N                  */
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "CrossP"
PetscReal CrossP(PetscReal M[3],PetscReal N[3], PetscInt m)
{
	if (m==0) return M[1]*N[2] - M[2]*N[1];
	if (m==1) return M[2]*N[0] - M[0]*N[2];
	if (m==2) return M[0]*N[1] - M[1]*N[0];
	return 0;
}

/* ------------------------------------------------------------------- */
/* Calculate the max of absolute values                                */
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "MaxAbs"
PetscReal MaxAbs(PetscReal a,PetscReal b)
{
	if (PetscAbsScalar(a) > PetscAbsScalar(b)) return PetscAbsScalar(a);
	else return PetscAbsScalar(b);
}

/* ------------------------------------------------------------------- */
/* Calculate the sum of absolute values                                */
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "SumAbs"
PetscReal SumAbs(PetscReal a,PetscReal b)
{
	return PetscAbsScalar(a)+PetscAbsScalar(b);
}

/* ------------------------------------------------------------------- */
/* Calculate the min of absolute values                                */
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "MinAbs"
PetscReal MinAbs(PetscReal a,PetscReal b)
{
	if (PetscAbsScalar(a) < PetscAbsScalar(b)) return PetscAbsScalar(a);
	else return PetscAbsScalar(b);
}

/* ------------------------------------------------------------------- */
/* Calculate the square of the norm 2 of P                             */
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "Norm2"
PetscReal Norm2(PetscReal M[3])
{
	return (M[0]*M[0]+M[1]*M[1]+M[2]*M[2]);
}

#undef __FUNCT__
#define __FUNCT__ "Average"
PetscErrorCode Average(Vec U,void* ctx)
{
	PetscErrorCode ierr;
	AppCtx         *user = (AppCtx*)ctx;
	DM             da = user->da;
	PetscInt       i,j,k,l,xs,ys,zs,xm,ym,zm;
	PetscReal      ****u,****u_;
	Vec            localU,localU_;
	
	PetscFunctionBegin;
	PetscInt rank;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMGetLocalVector(da,&localU );CHKERRQ(ierr);
	ierr = DMGetLocalVector(da,&localU_);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da,U,INSERT_VALUES,localU );CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(  da,U,INSERT_VALUES,localU );CHKERRQ(ierr);
	ierr = VecCopy(localU,localU_); CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(da,localU ,&u );CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(da,localU_,&u_);CHKERRQ(ierr);
	
	for (k=zs; k<zs+zm; k++)
		for (j=ys; j<ys+ym; j++)
			for (i=xs; i<xs+xm; i++)
				for (l=0; l<19; l++)
					u_[k][j][i][l] = (u[k+1][j][i][l] + u[k][j+1][i][l] + u[k][j][i+1][l] + u[k][j][i-1][l] + u[k][j-1][i][l] + u[k-1][j][i][l])/6.0;

	ierr = DMDAVecRestoreArrayDOF(da,localU ,&u );CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(da,localU_,&u_);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da,localU_,INSERT_VALUES,U);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(  da,localU_,INSERT_VALUES,U);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(da,&localU_);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(da,&localU);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "WeightedAverage"
PetscErrorCode WeightedAverage(Vec U,void* ctx)
{
	PetscErrorCode ierr;
	AppCtx         *user = (AppCtx*)ctx;
	DM             da = user->da;
	PetscInt       i,j,k,l,xs,ys,zs,xm,ym,zm;
	PetscReal      ****u,****u_;
	Vec            localU,localU_;
	
	PetscFunctionBegin;
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	ierr = DMGetLocalVector(da,&localU );CHKERRQ(ierr);
	ierr = DMGetLocalVector(da,&localU_);CHKERRQ(ierr);
	ierr = DMGlobalToLocalBegin(da,U,INSERT_VALUES,localU );CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(  da,U,INSERT_VALUES,localU );CHKERRQ(ierr);
	ierr = VecCopy(localU,localU_); CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(da,localU ,&u );CHKERRQ(ierr);
	ierr = DMDAVecGetArrayDOF(da,localU_,&u_);CHKERRQ(ierr);
	for (k=zs; k<zs+zm; k++)
		for (j=ys; j<ys+ym; j++)
			for (i=xs; i<xs+xm; i++)
				for (l=0; l<19; l++)
					u_[k][j][i][l] = (u[k+1][j][i][l] + u[k][j+1][i][l] + u[k][j][i+1][l] + 6*u[k][j][i][l] + u[k][j][i-1][l] + u[k][j-1][i][l] + u[k-1][j][i][l])/12.0;

	ierr = DMDAVecRestoreArrayDOF(da,localU ,&u );CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(da,localU_,&u_);CHKERRQ(ierr);
	ierr = DMLocalToGlobalBegin(da,localU_,INSERT_VALUES,U);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(  da,localU_,INSERT_VALUES,U);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(da,&localU_);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(da,&localU);CHKERRQ(ierr);
	
	PetscFunctionReturn(0);
}

