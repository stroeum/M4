/*
 * myCtx.c
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
  if (smoothing) ierr = WeightedAverage(F,user);
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

  // if (user->cnt == 0) exit(1);
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
  
  ierr = PetscGetTime(&user->t0);CHKERRQ(ierr);
  user->ldName = sprintf(user->dName,"ouput/");
  ierr = PetscOptionsGetString(PETSC_NULL,"-output_folder",user->dName,PETSC_MAX_PATH_LEN-1,PETSC_NULL);CHKERRQ(ierr);
  user->ldName = sprintf(user->vName,"viz_dir/");
  ierr = PetscOptionsGetString(PETSC_NULL,"-visualization_folder",user->vName,PETSC_MAX_PATH_LEN-1,PETSC_NULL);CHKERRQ(ierr);

  user->isInputFile = PETSC_FALSE;
  user->ldName = sprintf(user->InputFile,"X0.bin");
  ierr = PetscOptionsGetString(PETSC_NULL,"-input_file",user->InputFile,PETSC_MAX_PATH_LEN-1,&user->isInputFile);CHKERRQ(ierr);

  usrmnt->drawcontours = PETSC_FALSE;
  ierr = PetscOptionsHasName(PETSC_NULL,"-drawcontours",&usrmnt->drawcontours);CHKERRQ(ierr);
  user->maxsteps  = 100;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-ts_max_steps",&user->maxsteps,PETSC_NULL);CHKERRQ(ierr);
  user->bcType = 0;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-bc_type",&user->bcType,PETSC_NULL);CHKERRQ(ierr);
  
  ierr = PetscOptionsGetReal(PETSC_NULL,"-mi0",&user->mi[0],PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-mi1",&user->mi[1],PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-mi2",&user->mi[2],PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-me",&user->me,PETSC_NULL);CHKERRQ(ierr);
  
  user->un[0] = 0.0;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-WindNS",&user->un[0],PETSC_NULL);CHKERRQ(ierr);
  user->un[1] = 0.0;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-WindEW",&user->un[1],PETSC_NULL);CHKERRQ(ierr);
  user->un[2] = 0.0;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-WindUD",&user->un[2],PETSC_NULL);CHKERRQ(ierr);

  user->xtra_out = PETSC_FALSE;
  ierr = PetscOptionsHasName(PETSC_NULL,"-ExtraDiagnostics",&user->xtra_out);CHKERRQ(ierr);
  user->smoothing = PETSC_FALSE;
  ierr = PetscOptionsHasName(PETSC_NULL,"-smoothing",&user->smoothing);CHKERRQ(ierr);
  user->limiters = PETSC_FALSE;
  ierr = PetscOptionsHasName(PETSC_NULL,"-limiters",&user->limiters);CHKERRQ(ierr);
  user->blockers = PETSC_FALSE;
  ierr = PetscOptionsHasName(PETSC_NULL,"-blockers",&user->blockers);CHKERRQ(ierr);
  user->BfieldType = 0;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-B_field_type",&user->BfieldType,PETSC_NULL);CHKERRQ(ierr);
  user->B[0] = 0.0;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-B0",&user->B[0],PETSC_NULL);CHKERRQ(ierr);
  user->B[1] = 0.0;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-B1",&user->B[1],PETSC_NULL);CHKERRQ(ierr);
  user->B[2] = 0.0;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-B2",&user->B[2],PETSC_NULL);CHKERRQ(ierr);
  user->B[3] = 0.0;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-B3",&user->B[3],PETSC_NULL);CHKERRQ(ierr);
  
  ierr = PetscOptionsGetReal(PETSC_NULL,"-Rmars",&user->rM,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-Mmars",&user->mM,PETSC_NULL);CHKERRQ(ierr);
  
  user->inXmin = 0.0;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-in_x_min",&user->inXmin,PETSC_NULL);CHKERRQ(ierr);
  user->inXmax = 1.0;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-in_x_max",&user->inXmax,PETSC_NULL);CHKERRQ(ierr);
  user->inYmin = 0.0;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-in_y_min",&user->inYmin,PETSC_NULL);CHKERRQ(ierr);
  user->inYmax = 1.0;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-in_y_max",&user->inYmax,PETSC_NULL);CHKERRQ(ierr);
  user->inZmin = 0.0;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-in_z_min",&user->inZmin,PETSC_NULL);CHKERRQ(ierr);
  if(user->inZmin<0) user->inZmin = 0.0;
  user->inZmax = 1.0;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-in_z_max",&user->inZmax,PETSC_NULL);CHKERRQ(ierr);
  user->outXmin = user->inXmin;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-out_x_min",&user->outXmin,PETSC_NULL);CHKERRQ(ierr);
  user->outXmax = user->inXmax;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-out_x_max",&user->outXmax,PETSC_NULL);CHKERRQ(ierr);
  user->outYmin = user->inYmin;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-out_y_min",&user->outYmin,PETSC_NULL);CHKERRQ(ierr);
  user->outYmax = user->inYmax;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-out_y_max",&user->outYmax,PETSC_NULL);CHKERRQ(ierr);
  user->outZmin = user->inZmin;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-out_z_min",&user->outZmin,PETSC_NULL);CHKERRQ(ierr);
  if(user->outZmin<0) user->outZmin = 0.0;
  user->outZmax = user->inZmax;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-out_z_max",&user->outZmax,PETSC_NULL);CHKERRQ(ierr);

  if(user->outXmin >= user->inXmin) user->outXmin = user->inXmin;
  if(user->outXmax <= user->inXmax) user->outXmax = user->inXmax;
  if(user->outYmin >= user->inYmin) user->outYmin = user->inYmin;
  if(user->outYmax <= user->inYmax) user->outYmax = user->inYmax;
  if(user->outZmin >= user->inZmin) user->outZmin = user->inZmin;
  if(user->outZmax <= user->inZmax) user->outZmax = user->inZmax;
  
  ierr = PetscOptionsGetReal(PETSC_NULL,"-vizbox_x_min",&user->vizbox[0],PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-vizbox_x_max",&user->vizbox[1],PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-vizbox_y_min",&user->vizbox[2],PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-vizbox_y_max",&user->vizbox[3],PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-vizbox_z_min",&user->vizbox[4],PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-vizbox_z_max",&user->vizbox[5],PETSC_NULL);CHKERRQ(ierr);
  
  user->viz_dstep = 0;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-viz_dstep",&user->viz_dstep,PETSC_NULL);CHKERRQ(ierr);
  
  mx_in = 2;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-da_grid_x",&mx_in,PETSC_NULL);CHKERRQ(ierr);
  my_in = 2;
  ierr = PetscOptionsGetInt(PETSC_NULL,"-da_grid_y",&my_in,PETSC_NULL);CHKERRQ(ierr);
  mz_in = 2; 
  ierr = PetscOptionsGetInt(PETSC_NULL,"-da_grid_z",&mz_in,PETSC_NULL);CHKERRQ(ierr);
  
  user->ti = 0.0; 
  ierr = PetscOptionsGetReal(PETSC_NULL,"-ts_init_time",&user->ti,PETSC_NULL);CHKERRQ(ierr);
  user->tf = 1.0;
  ierr = PetscOptionsGetReal(PETSC_NULL,"-ts_max_time",&user->tf,PETSC_NULL);CHKERRQ(ierr);

  
  if(user->viz_dstep<=0 || (PetscInt)(user->viz_dstep)!=user->viz_dstep){
    PetscPrintf(PETSC_COMM_WORLD,"ERROR: Invalid value of viz_dstep: viz_dstep=1");
    user->viz_dstep = 1;
  }
  
  user->gM = G*user->mM/(user->rM*user->rM);
  
  user->Lx = user->inXmax-user->inXmin;
  user->Ly = user->inYmax-user->inYmin;
  user->Lz = user->inZmax-user->inZmin;
  user->L  = 1.0e3;//user->Lz;

  dx = user->Lx/(PetscReal)(mx_in-1);
  dy = user->Ly/(PetscReal)(my_in-1);
  dz = user->Lz/(PetscReal)(mz_in-1);

  mx_lo = (PetscInt)ceil(log((pct-1.0)*(user->inXmin-user->outXmin)/(pct*dx) + 1)/log(pct)); 
  my_lo = (PetscInt)ceil(log((pct-1.0)*(user->inYmin-user->outYmin)/(pct*dy) + 1)/log(pct));
  mz_lo = (PetscInt)ceil(log((pct-1.0)*(user->inZmin-user->outZmin)/(pct*dz) + 1)/log(pct));

  mx_hi = (PetscInt)floor(log((pct-1.0)*(user->outXmax-user->inXmax)/(pct*dx) + 1)/log(pct)); 
  my_hi = (PetscInt)floor(log((pct-1.0)*(user->outYmax-user->inYmax)/(pct*dy) + 1)/log(pct));
  mz_hi = (PetscInt)floor(log((pct-1.0)*(user->outZmax-user->inZmax)/(pct*dz) + 1)/log(pct));

  // Calculate dimension of the total system from outerbound to outerbound //
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
    //PetscPrintf(PETSC_COMM_WORLD,"X = %12.3f\n",user->outXmin + user->x[i]*user->L);
  }
  user->outYmin = user->inYmin - pct*dy * (pow(pct,my_lo)-1) / (pct-1);
  for (j=0; j<user->my; j++) {
    if (j < my_lo) Y = user->inYmin - pct*dy * (pow(pct,my_lo-j)-1) / (pct-1);
    else if (j >= my_lo && j < my_lo+my_in) Y = user->inYmin + (j-my_lo)*dy;
    else Y = user->inYmax + pct*dy * (pow(pct,j-my_lo-my_in+1)-1) / (pct-1);
    user->y[j] = (Y - user->outYmin)/user->L;
    //PetscPrintf(PETSC_COMM_WORLD,"Y = %12.3f\n",user->outYmin + user->y[j]*user->L);
  }
  user->outZmin = user->inZmin - pct*dz * (pow(pct,mz_lo)-1) / (pct-1);
  for (k=0; k<user->mz; k++) {
    if (k < mz_lo) Z = user->inZmin - pct*dz * (pow(pct,mz_lo-k)-1) / (pct-1);
    else if (k >= mz_lo && k < mz_lo+mz_in) Z = user->inZmin + (k-mz_lo)*dz;
    else Z = user->inZmax + pct*dz * (pow(pct,k-mz_lo-mz_in+1)-1) / (pct-1);
    user->z[k] = (Z - user->outZmin)/user->L;
    //PetscPrintf(PETSC_COMM_WORLD,"Z = %12.3f\n",user->outZmin + user->z[k]*user->L);
  }

  user->tau = PetscSqrtScalar(user->L/user->gM);
  user->n0  = user->me/(qe*qe*user->L*user->L*mu0);
  user->v0  = PetscSqrtScalar(user->gM*user->L);
  user->p0  = (user->me*user->me*user->gM)/(qe*qe*user->L*mu0);
  user->B0  = user->me/(qe*user->tau);
  
  /* dt = min(dt,CFL) */
  user->dt    = 1.0e-9/user->tau;
  user->eps.n = 1.0e3/user->n0; // min density allowed in the domain: 10^-3 cm-3 = 1e3 m^-3
  user->eps.p = 1e-15/user->p0; // min pressure allowed in the domain
  user->eps.v = .01*299742458/user->v0; // MAX velocity allowed in the domain: 1% of the speed of light

  // Retrieve information from initial conditions //
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
    if(pdot)
    {
      size_t len = pdot - InputFile -1; // string length
      pstep = malloc(len);
      memcpy(pstep, &InputFile[1],len);
      pstep[len] = '\0';
      user->istep = atoi(pstep);
      PetscPrintf(PETSC_COMM_WORLD,"Initial timestep %D\n",user->istep);

      /*
       * Retrieve initial time from step number and file t.dat containing all the time sequence.
       * 1. Retrieve the line number
       * 2. Scan the file t.dat until this line
       * 3. Store the time in the variable t
       * 4. Store the initial time in user->ti
       */
      PetscInt  tg_line;
      PetscInt  line=1;
      tg_line = 1 + user->istep/user->viz_dstep;
      PetscPrintf(PETSC_COMM_WORLD,"target line number %d\n",tg_line);
      FILE *pFile;
      float t; 
      pFile = fopen("output/t.dat","r");
      while(pFile && line<=tg_line) {
        fscanf(pFile,"%e",&t);
        //PetscPrintf(PETSC_COMM_WORLD,"line %d\tti = %e\n",line, t);
        line++;
      }
      fclose(pFile);
      //user->ti = user->istep*user->dt;
    } else {
      PetscPrintf(PETSC_COMM_WORLD,"Unreadable input file\n");
      exit(1);
    }
  }

  if(user->vizbox[0]<user->outXmin || user->vizbox[1]>user->outXmax || user->vizbox[2]<user->outYmin || user->vizbox[3]>user->outYmax || user->vizbox[4]<user->outZmin || user->vizbox[5]>user->outZmax){
    PetscPrintf(PETSC_COMM_WORLD,"WARNING: Visualization area greater than simulation domain:\n---> the vizualization box is adjusted to fit the maximum dimension of the simulation domain.\n");
    if (user->vizbox[0]<user->outXmin) user->vizbox[0] = user->outXmin;
    if (user->vizbox[1]>user->outXmax) user->vizbox[1] = user->outXmax;
    if (user->vizbox[2]<user->outYmin) user->vizbox[2] = user->outYmin;
    if (user->vizbox[3]>user->outYmax) user->vizbox[3] = user->outYmax;
    if (user->vizbox[4]<user->outZmin) user->vizbox[4] = user->outZmin;
    if (user->vizbox[5]>user->outZmax) user->vizbox[5] = user->outZmax;
  }

  // Read input files //
  if(rank==0) {
    ReadTable(user->RefProf,Np_REF,"input/Profiles.dat");
    ReadTable(user->RefPart,4,"input/Partition.dat");
  }
  MPI_Bcast(&(user->RefProf[0][0]),Np_REF*Nz_REF, MPIU_REAL, 0, PETSC_COMM_WORLD);
  MPI_Bcast(&(user->RefPart[0][0]),     4*Nz_REF, MPIU_REAL, 0, PETSC_COMM_WORLD);

  // Display diagnostics //
  PetscPrintf(PETSC_COMM_WORLD,"mx = [%i %i %i]:%i\nmy = [%i %i %i]:%i\nmz = [%i %i %i]:%i\n",mx_lo,mx_in,mx_hi,user->mx, my_lo,my_in,my_hi,user->my, mz_lo,mz_in,mz_hi,user->mz); 
  PetscPrintf(PETSC_COMM_WORLD,"      O2+\t\tCO2+\t\tO+\t\te\n");
  PetscPrintf(PETSC_COMM_WORLD,"m   = [%2.3e\t%2.3e\t%2.3e\t%2.3e]\n",user->mi[0],user->mi[1],user->mi[2],user->me);
  PetscPrintf(PETSC_COMM_WORLD,"vo  = [%2.3e]\n",user->v0);
  PetscPrintf(PETSC_COMM_WORLD,"tau = [%2.3e]\n",user->tau);
  PetscPrintf(PETSC_COMM_WORLD,"n0  = [%2.3e]\n",user->n0);
  PetscPrintf(PETSC_COMM_WORLD,"p0  = [%2.3e]\n",user->p0);
  PetscPrintf(PETSC_COMM_WORLD,"B0  = [%2.3e]\n",user->B0);
  PetscPrintf(PETSC_COMM_WORLD,"E0  = [%2.3e]\n",user->v0*user->B0);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "OutputData"
PetscErrorCode OutputData(void* ptr)
{
  PetscErrorCode ierr;
  AppCtx         *user = (AppCtx*)ptr;
  PetscInt       istep=user->istep,maxsteps = user->maxsteps,vizdstep = user->viz_dstep;
  PetscBool      xtra_out = user->xtra_out;
  PetscReal      Xmin = user->outXmin, Ymin = user->outYmin, Zmin = user->outZmin;
  PetscReal      n0 = user->n0, v0 = user->v0, p0 = user->p0, B0 = user->B0;
  PetscReal      *x = user->x, *y = user->y, *z = user->z;
  PetscInt       mx = user->mx, my = user->my, mz = user->mz;
  DM             da = user->da;
  DM             db = user->db;
  PetscReal      L = user->L;
  dvi            d = user->d;
  svi            s = user->s;
  PetscReal      vizbox[6] = {
    user->vizbox[0],user->vizbox[1],
    user->vizbox[2],user->vizbox[3],
    user->vizbox[4],user->vizbox[5]
  };

  Vec            U,V;
  PetscReal      ****u,****v;
  FILE           *pFile;
  char           fName[PETSC_MAX_PATH_LEN];
  PetscViewer    fViewer;
  PetscInt       rank,step,flag;
  PetscInt       i,j,k,l,m,xs,ys,zs,xm,ym,zm;
  PetscInt       imin,imax,jmin,jmax,kmin,kmax; 
  PetscReal      ne,ni[3],ve[3],vi[3][3],pe,pi[3],J[3],E[3],B[3]; 
  PetscReal      X,Y,Z;

  PetscFunctionBegin;

  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  PetscPrintf(PETSC_COMM_WORLD,"box: [%2.1f,%2.1f,%2.1f] to [%2.1f,%2.1f,%2.1f] (km)\n",vizbox[0]/1e3,vizbox[2]/1e3,vizbox[4]/1e3,vizbox[1]/1e3,vizbox[3]/1e3,vizbox[5]/1e3);
  
  // Get local grid boundaries //
  ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  for (step=istep; step<=maxsteps; step++) {

    if(step%vizdstep==0) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"timestep %D\n",step);CHKERRQ(ierr);

      // Read binary file containing the solution //
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
      
      // Get pointers to vector data //
      ierr = DMDAVecGetArrayDOF(da,U,&u);CHKERRQ(ierr);
      if(xtra_out) {ierr = DMDAVecGetArrayDOF(db,V,&v);CHKERRQ(ierr);}
      
      // Store the solution limited to the visualization box into an ascii file //
      sprintf(fName,"%s/X%d.dat",user->vName,step);
      flag  = access(fName,W_OK);
      if (flag==0) pFile = fopen(fName,"a");
      else         pFile = fopen(fName,"w");
      
      imin = mx;
      imax = 0;
      jmin = my;
      jmax = 0;
      kmin = mz;
      kmax = 0;

      for (i=xs; i<xs+xm; i++) {
        X = Xmin + x[i]*L;
        if((float)X>=vizbox[0] && (float)X<=vizbox[1]){
          for (j=ys; j<ys+ym; j++) {
            Y = Ymin + y[j]*L;
            if((float)Y>=vizbox[2] && (float)Y<=vizbox[3]){
              for (k=zs; k<zs+zm; k++) {
                Z = Zmin + z[k]*L;
                if((float)Z>=vizbox[4] && (float)Z<=vizbox[5]){
                  if (i<imin) imin = i;
                  if (i>imax) imax = i;
                  if (j<jmin) jmin = j;
                  if (j>jmax) jmax = j;
                  if (k<kmin) kmin = k;
                  if (k>kmax) kmax = k;
                  
                  if(xtra_out) {
                    ne = (u[k][j][i][d.ni[O2p]]+u[k][j][i][d.ni[CO2p]]+u[k][j][i][d.ni[Op]])*n0;
                    ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"%2.3e\t",ne);CHKERRQ(ierr);
                  }
                  for (l=0; l<3; l++) {
                    ni[l] = u[k][j][i][d.ni[l]]*n0;
                    ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"%2.3e\t",ni[l]);CHKERRQ(ierr);
                  }
                  if(xtra_out) {
                    for (m=0; m<3; m++) {
                      ve[m] = v[k][j][i][s.ve[m]]*v0;
                      ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"%2.3e\t",ve[m]);CHKERRQ(ierr);
                    }
                  }
                  for (l=0; l<3; l++) {
                    for (m=0; m<3; m++) {
                      vi[l][m] = u[k][j][i][d.vi[l][m]]*v0;
                      ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"%2.3e\t",vi[l][m]);CHKERRQ(ierr);
                    }
                  }
                  pe = u[k][j][i][d.pe]*p0;
                  ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"%2.3e\t",pe);CHKERRQ(ierr);
                  for (l=0; l<3; l++) {
                    pi[l] = u[k][j][i][d.pi[l]]*p0;
                    ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"%2.3e\t",pi[l]);CHKERRQ(ierr);
                  }
                  if(xtra_out) {
                    for (m=0; m<3; m++){
                      J[m] = v[k][j][i][s.J[m]]*qe*n0*v0;
                      ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"%2.3e\t",J[m]);CHKERRQ(ierr);
                    }
                    for (m=0; m<3; m++){
                      E[m] = v[k][j][i][s.E[m]]*v0*B0;
                      ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"%2.3e\t",E[m]);CHKERRQ(ierr);
                    }
                  }
                  for (m=0; m<3; m++){
                    B[m] = u[k][j][i][d.B[m]]*B0;
                    ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"%2.3e\t",B[m]);CHKERRQ(ierr);
                  }
                  ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"\n");CHKERRQ(ierr);
                }
              }
            }
          }
        }
      }
      fclose(pFile);

      m = MPI_Allreduce(MPI_IN_PLACE,&imin,1,MPIU_INT,MPIU_MIN,PETSC_COMM_WORLD);
      m = MPI_Allreduce(MPI_IN_PLACE,&imax,1,MPIU_INT,MPIU_MAX,PETSC_COMM_WORLD);
      m = MPI_Allreduce(MPI_IN_PLACE,&jmin,1,MPIU_INT,MPIU_MIN,PETSC_COMM_WORLD);
      m = MPI_Allreduce(MPI_IN_PLACE,&jmax,1,MPIU_INT,MPIU_MAX,PETSC_COMM_WORLD);
      m = MPI_Allreduce(MPI_IN_PLACE,&kmin,1,MPIU_INT,MPIU_MIN,PETSC_COMM_WORLD);
      m = MPI_Allreduce(MPI_IN_PLACE,&kmax,1,MPIU_INT,MPIU_MAX,PETSC_COMM_WORLD);

      sprintf(fName,"%s/X%d.general",user->vName,step);
      flag  = access(fName,W_OK);
      if (flag==0) pFile = fopen(fName,"a");
      else         pFile = fopen(fName,"w");

      if(!xtra_out){
        ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"file = X%d.dat\n",step);CHKERRQ(ierr);
        ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"grid = %d %d %d\n",imax-imin+1,jmax-jmin+1,kmax-kmin+1);CHKERRQ(ierr);
        ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"format = ascii\n");CHKERRQ(ierr);
        ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"interleaving = field\n");CHKERRQ(ierr);
        ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"field = nO2p,nCO2p,nOp,vO2p,vCO2p,vOp,pO2p,pCO2p,pOp,pe,B\n");CHKERRQ(ierr);
        ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"type = float float float float float float float float float float float\n");CHKERRQ(ierr);
        ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"structure = scalar scalar scalar 3-vector 3-vector 3-vector scalar scalar scalar scalar 3-vector\n");CHKERRQ(ierr);
        ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"dependency = positions,positions,positions,positions,positions,positions,positions,positions,positions,positions,positions\n");CHKERRQ(ierr);
        //ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"positions = regular, regular, regular, %2.3e, %2.3e, %2.3e, %2.3e, %2.3e, %2.3e\n", (Xmin+x[imin]*L)*1e-3,(x[imax]-x[imin])/(imax-imin)*L*1e-3,(Ymin+y[jmin]*L)*1e-3,(y[jmax]-y[jmin])/(jmax-jmin)*L*1e-3,(Zmin+z[kmin]*L)*1e-3,(z[kmax]-z[kmin])/(kmax-kmin)*L*1e-3);CHKERRQ(ierr);
        ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"positions = irregular, irregular, irregular");CHKERRQ(ierr);
        for (i=imin; i<=imax; i++)  ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,", %2.3e",(Xmin+x[i]*L)*1e-3);CHKERRQ(ierr);
        for (j=jmin; j<=jmax; j++)  ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,", %2.3e",(Ymin+y[j]*L)*1e-3);CHKERRQ(ierr);
        for (k=kmin; k<=kmax; k++)  ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,", %2.3e",(Zmin+z[k]*L)*1e-3);CHKERRQ(ierr);
        ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"\n\nend");CHKERRQ(ierr);
      } else {
        ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"file = X%d.dat\n",step);CHKERRQ(ierr);
        ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"grid = %d %d %d\n",imax-imin+1,jmax-jmin+1,kmax-kmin+1);CHKERRQ(ierr);
        ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"format = ascii\n");CHKERRQ(ierr);
        ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"interleaving = field\n");CHKERRQ(ierr);
        ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"field = ne,nO2p,nCO2p,nOp,ve,vO2p,vCO2p,vOp,pe,pO2p,pCO2p,pOp,J,E,B\n");CHKERRQ(ierr);
        ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"type = float float float float float float float float float float float float float float float\n");CHKERRQ(ierr);
        ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"structure = scalar scalar scalar scalar 3-vector 3-vector 3-vector 3-vector scalar scalar scalar scalar 3-vector 3-vector 3-vector\n");CHKERRQ(ierr);
        ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"dependency = positions,positions,positions,positions,positions,positions,positions,positions,positions,positions,positions,positions,positions,positions,positions\n");CHKERRQ(ierr);
        //ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"positions = regular, regular, regular, %2.3e, %2.3e, %2.3e, %2.3e, %2.3e, %2.3e\n", (Xmin+x[imin]*L)*1e-3,(x[imax]-x[imin])/(imax-imin)*L*1e-3,(Ymin+y[jmin]*L)*1e-3,(y[jmax]-y[jmin])/(jmax-jmin)*L*1e-3,(Zmin+z[kmin]*L)*1e-3,(z[kmax]-z[kmin])/(kmax-kmin)*L*1e-3);CHKERRQ(ierr);
        ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"positions = irregular, irregular, irregular");CHKERRQ(ierr);
        for (i=imin; i<=imax; i++)  ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,", %2.3e",(Xmin+x[i]*L)*1e-3);CHKERRQ(ierr);
        for (j=jmin; j<=jmax; j++)  ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,", %2.3e",(Ymin+y[j]*L)*1e-3);CHKERRQ(ierr);
        for (k=kmin; k<=kmax; k++)  ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,", %2.3e",(Zmin+z[k]*L)*1e-3);CHKERRQ(ierr);
        ierr = PetscFPrintf(PETSC_COMM_WORLD,pFile,"\n\nend");CHKERRQ(ierr);
      }
      fclose(pFile);
      ierr = DMDAVecRestoreArrayDOF(da,U,&u);CHKERRQ(ierr);
      ierr = VecDestroy(&U);CHKERRQ(ierr);
      if(xtra_out) {
        ierr = DMDAVecRestoreArrayDOF(db,V,&v);CHKERRQ(ierr);
        ierr = VecDestroy(&V);CHKERRQ(ierr);
      }
    }
  }
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
  PetscReal      old_dt,new_dt;

  PetscFunctionBegin;
  ierr = TSGetApplicationContext(ts,(void *)&user);CHKERRQ(ierr);
  ierr = TSGetTimeStep(ts,&old_dt);CHKERRQ(ierr);
  new_dt = user->dt;
  //PetscPrintf(PETSC_COMM_WORLD,"dt = %3.5e %3.5e\n", old_dt*user->tau,new_dt*user->tau);
  //ierr = TSSetTimeStep(ts,PetscMin(old_dt,new_dt));CHKERRQ(ierr);
  //ierr = TSSetTimeStep(ts,new_dt);CHKERRQ(ierr);
  ierr = TSSetTimeStep(ts,user->dt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
/* Calculate the magnetic field of multiple vertical dipoles           */
/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "ARCADES"
/*
 * Arcades
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
        s = k*Ns*Ms + j*Ms + i;
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

  /*
  if(rank==0) {
  for (k=-1;k<=1;k++) for (j=-1;j<=1;j++) for (i=-1;i<=1;i++) 
    if (! ( (i==-1 && j==-1) || (i==-1 && k==-1) || (j==-1 && k==-1) ) ) {
      PetscPrintf(PETSC_COMM_SELF,"@[%2d][%2d][%2d]  nO2+=%2.5e;  nCO2+=%2.5e;  nO+=%2.5e;  vO2+=[%2.5e %2.5e %2.5e];  vCO2+=[%2.5e %2.5e %2.5e];  vO+=[%2.5e %2.5e %2.5e];  pO2+=%2.5e;  pCO2+=%2.5e;  pO+=%2.5e;  pe=%2.5e;  B=[%2.5e %2.5e %2.5e]\n",i,j,k,u[k][j][i][0]*user->n0,u[k][j][i][1]*user->n0,u[k][j][i][2]*user->n0,u[k][j][i][3]*user->v0,u[k][j][i][4]*user->v0,u[k][j][i][5]*user->v0,u[k][j][i][6]*user->v0,u[k][j][i][7]*user->v0,u[k][j][i][8]*user->v0,u[k][j][i][9]*user->v0,u[k][j][i][10]*user->v0,u[k][j][i][11]*user->v0,u[k][j][i][12]*user->p0,u[k][j][i][13]*user->p0,u[k][j][i][14]*user->p0,u[k][j][i][15]*user->p0,u[k][j][i][16]*user->B0,u[k][j][i][17]*user->B0,u[k][j][i][18]*user->B0);
    }
  }
  if (rank == 0) { 
  for (k=-1;k<=1;k++) for (j=-1;j<=1;j++) for (i=-1;i<=1;i++) 
    if (! ( (i==-1 && j==-1) || (i==-1 && k==-1) || (j==-1 && k==-1) ) ) {
      PetscPrintf(PETSC_COMM_SELF,"@[%2d][%2d][%2d] _nO2+=%2.5e; _nCO2+=%2.5e; _nO+=%2.5e; _vO2+=[%2.5e %2.5e %2.5e]; _vCO2+=[%2.5e %2.5e %2.5e]; _vO+=[%2.5e %2.5e %2.5e]; _pO2+=%2.5e; _pCO2+=%2.5e; _pO+=%2.5e; _pe=%2.5e; _B=[%2.5e %2.5e %2.5e]\n",i,j,k,u_[k][j][i][0]*user->n0,u_[k][j][i][1]*user->n0,u_[k][j][i][2]*user->n0,u_[k][j][i][3]*user->v0,u_[k][j][i][4]*user->v0,u_[k][j][i][5]*user->v0,u_[k][j][i][6]*user->v0,u_[k][j][i][7]*user->v0,u_[k][j][i][8]*user->v0,u_[k][j][i][9]*user->v0,u_[k][j][i][10]*user->v0,u_[k][j][i][11]*user->v0,u_[k][j][i][12]*user->p0,u_[k][j][i][13]*user->p0,u_[k][j][i][14]*user->p0,u_[k][j][i][15]*user->p0,u_[k][j][i][16]*user->B0,u_[k][j][i][17]*user->B0,u_[k][j][i][18]*user->B0);
    }
  }
  if (rank == 0) { 
    i=0;
    j=0;
    k=0;
    PetscPrintf(PETSC_COMM_SELF,"@[%2d][%2d][%2d] avnO2p=%2.5e; avnCO2p=%2.5e; avnOp=%2.5e; avvO2p=[%2.5e %2.5e %2.5e]; avvCO2p=[%2.5e %2.5e %2.5e]; avvOp=[%2.5e %2.5e %2.5e]; avpO2p=%2.5e; avpCO2p=%2.5e; avpOp=%2.5e; avpe=%2.5e; avB=[%2.5e %2.5e %2.5e]\n",i,j,k,u_[k][j][i][0]*user->n0,u_[k][j][i][1]*user->n0,u_[k][j][i][2]*user->n0,u_[k][j][i][3]*user->v0,u_[k][j][i][4]*user->v0,u_[k][j][i][5]*user->v0,u_[k][j][i][6]*user->v0,u_[k][j][i][7]*user->v0,u_[k][j][i][8]*user->v0,u_[k][j][i][9]*user->v0,u_[k][j][i][10]*user->v0,u_[k][j][i][11]*user->v0,u_[k][j][i][12]*user->p0,u_[k][j][i][13]*user->p0,u_[k][j][i][14]*user->p0,u_[k][j][i][15]*user->p0,u_[k][j][i][16]*user->B0,u_[k][j][i][17]*user->B0,u_[k][j][i][18]*user->B0);
  }
  */

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

