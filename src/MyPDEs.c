#include "MyPDEs.h"

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormInitialSolution"
/* 
 FormInitiatialSolution - Initial conditions for u.
 
 Input Parameters:
 .  ts - the TS context
 .  U - input vector
 .  ctx - optional user-defined context, as set by SNESSetFunction()
 
 Output Parameter:
 .  F - function vector
 */
PetscErrorCode FormInitialSolution(Vec U, void* ctx)
{
	PetscErrorCode ierr;
	AppCtx         *user = (AppCtx*)ctx;
	tolerance      eps = user->eps;
	DM             da = user->da;
	dvi            d = user->d;
	PetscInt       Btype = user->BfieldType;
	PetscReal      *x=user->x, *y=user->y, *z=user->z;
	PetscReal      L=user->L,Lz;
	PetscReal      Xmin=user->outXmin, Ymin=user->outYmin, Zmin=user->outZmin;
	PetscReal      inZmax=user->inZmax;
	PetscReal      n0=user->n0, v0=user->v0, p0=user->p0, B0=user->B0;
	PetscReal      Bo=user->B[0], a=user->B[1], b=user->B[2], c=user->B[3];
	PetscReal      ui[3]={user->ui[0],user->ui[1],user->ui[2]}; // neutral wind
	
	PetscBool      vDamping = user->vDamping;
	PetscReal      lambda=user->lambda;
	PetscReal      XL=user->inXmin,XU=user->inXmax;
	PetscReal      YL=user->inYmin,YU=user->inYmax;
	PetscReal      ZL=user->ZL    ,ZU=user->ZU    ;
	
	PetscInt       i,j,k,l,m,xs,ys,zs,xm,ym,zm;
	PetscReal      ****u;
	PetscReal      X,Y,Z;
	PetscReal      Te,Ti;
	PetscReal      nio[3],neo,pio[3],peo;
	PetscInt       rank;
	PetscViewer    fViewer;
	char           fName[PETSC_MAX_PATH_LEN];
	PetscBool      flag;
	
	PetscFunctionBegin;
	PetscPrintf(PETSC_COMM_WORLD,"\nForm Initial Conditions: ...\n");
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	Lz = user->outZmax-user->inZmin;
	if (user->isInputFile) {
		// Read binary file containing the solution
		sprintf(fName,"%s/%s",user->dName,user->InputFile);
		PetscPrintf(PETSC_COMM_WORLD,"Init from %s: ...\n",fName);
		flag = access(fName,F_OK);
		if (flag != 0) {
			PetscPrintf(PETSC_COMM_WORLD,"Unreadable input file\n");
			exit(1);
		}
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,fName,FILE_MODE_READ, &fViewer); CHKERRQ(ierr);
		ierr = VecLoad(U,fViewer); CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&fViewer);CHKERRQ(ierr);   
		PetscPrintf(PETSC_COMM_WORLD,"Init from %s: DONE\n",fName);
	} else { 
		// Map U
		ierr = DMDAVecGetArrayDOF(da,U,&u);CHKERRQ(ierr);
		
		// Get local grid boundaries
		ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
		
		for (k=zs; k<zs+zm; k++) {
			Z = Zmin + z[k]*L;
			for (j=ys; j<ys+ym; j++) {
				Y = Ymin + y[j]*L;
				for (i=xs; i<xs+xm; i++) {
					X = Xmin + x[i]*L;
					
					Te = Interpolate(user->RefProf, 8 ,Z, lin_flat);
					Ti = Interpolate(user->RefProf, 9 ,Z, lin_flat);
					neo = Interpolate(user->RefProf, 7 ,Z, lin_exp);
					assert(Te>0);
					assert(Ti>0);
					assert(neo>0);

					nio[O2p]  = Interpolate(user->RefPart, O2p ,Z, lin_flat)*neo;
					nio[CO2p] = Interpolate(user->RefPart, CO2p,Z, lin_flat)*neo;
					nio[Op]   = Interpolate(user->RefPart, Op  ,Z, lin_flat)*neo;

					if(nio[Op] <10) nio[Op] = 10;

					for (l=0; l<3; l++) {
						if (!(nio[l] > 0)) nio[l] = eps.n*n0;
                        
						pio[l] = nio[l]*kB*Ti;
						assert(pio[l]>0);
					}

					peo = neo*kB*Te;
					assert(peo>0);
					if(Btype==0) {
						u[k][j][i][d.B[0]] = Interpolate(user->RefProf, 0 ,Z, lin_flat)/B0;
						u[k][j][i][d.B[1]] = Interpolate(user->RefProf, 1 ,Z, lin_flat)/B0;
						u[k][j][i][d.B[2]] = Interpolate(user->RefProf, 2 ,Z, lin_flat)/B0;
					} else if (Btype==1) {
						u[k][j][i][d.B[0]] = 0.0/B0;
						u[k][j][i][d.B[1]] = 0.0/B0;
						u[k][j][i][d.B[2]] = Bo /B0;
					} else if (Btype==2) {
						u[k][j][i][d.B[0]] = 0.0/B0;
						u[k][j][i][d.B[1]] = Bo /B0;
						u[k][j][i][d.B[2]] = 0.0/B0;
					} else if (Btype==3) {
						u[k][j][i][d.B[0]] = V_Dipole(Bo,0,0,a,X,Y,Z,0)/B0;
						u[k][j][i][d.B[1]] = V_Dipole(Bo,0,0,a,X,Y,Z,1)/B0;
						u[k][j][i][d.B[2]] = V_Dipole(Bo,0,0,a,X,Y,Z,2)/B0;
					} else if (Btype==4) {
						u[k][j][i][d.B[0]] = H_Dipole(Bo,0,0,a,X,Y,Z,0)/B0;
						u[k][j][i][d.B[1]] = H_Dipole(Bo,0,0,a,X,Y,Z,1)/B0;
						u[k][j][i][d.B[2]] = H_Dipole(Bo,0,0,a,X,Y,Z,2)/B0;
					} else if (Btype==5) {
						if (Z<=inZmax) {
							u[k][j][i][d.B[0]] = 0.0/B0;
							u[k][j][i][d.B[1]] = 0.0/B0;
							u[k][j][i][d.B[2]] = Bo /B0;
						} else {
							u[k][j][i][d.B[0]] = 0.0/B0;
							u[k][j][i][d.B[1]] = 0.0/B0;
							u[k][j][i][d.B[2]] = Bo /B0 * (1-erf((Z-3*a-inZmax)/a))/2;
						}
					} else if (Btype==6) {
						if (Z<=inZmax) {
							u[k][j][i][d.B[0]] = 0.0/B0;
							u[k][j][i][d.B[1]] = 0.0/B0;
							u[k][j][i][d.B[2]] = Bo /B0;
						} else {
							u[k][j][i][d.B[0]] = 0.0/B0;
							u[k][j][i][d.B[1]] = 0.0/B0;
							u[k][j][i][d.B[2]] = Bo /B0 * pow(inZmax/Z,3);
						}
					} else if (Btype==7) {
						u[k][j][i][d.B[0]] = Arcades(X,Y,Z,0)/B0;
						u[k][j][i][d.B[1]] = Arcades(X,Y,Z,1)/B0;
						u[k][j][i][d.B[2]] = Arcades(X,Y,Z,2)/B0;
					} else if (Btype==8) {
						u[k][j][i][d.B[0]] = MultiArcades(X,Y,Z,0)/B0;
						u[k][j][i][d.B[1]] = MultiArcades(X,Y,Z,1)/B0;
						u[k][j][i][d.B[2]] = MultiArcades(X,Y,Z,2)/B0;
					}
					
					for (l=0; l<3; l++) {
						u[k][j][i][d.ni[l]] = nio[l]/n0;
						for (m=0; m<3; m++) {
							if (vDamping) { //((1.0-tanh(Z-Lz)/(Lz/12.0))/2.0);
								u[k][j][i][d.vi[l][m]] = ui[m]/v0 * 
								.5*(1+erf((X-XL)/lambda)) * .5*(1-erf((X-XU)/lambda)) * 
								.5*(1+erf((Y-YL)/lambda)) * .5*(1-erf((Y-YU)/lambda)) * 
								.5*(1+erf((Z-ZL)/lambda)) * .5*(1-erf((Z-ZU)/lambda)) ;
							} else {
								u[k][j][i][d.vi[l][m]] = ui[m]/v0;
							}
						}
						u[k][j][i][d.pi[l]] = pio[l]/p0;
					}
					u[k][j][i][d.pe] = peo/p0;
				}
			}
		}
		
		ierr = DMDAVecRestoreArrayDOF(da,U,&u);CHKERRQ(ierr);
	}
	
	PetscPrintf(PETSC_COMM_WORLD,"Form Initial Conditions: DONE\n");
	PetscFunctionReturn(0); 
} 

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormBCu"
/* 
 FormBoundaryConditions - Evaluates nonlinear function @ BC (outside the domain), F(u).
 * T: Top, B: Bottom, N: North, S: South, E: East, W: West
 * 8 corners: NWD, SWD, NED, SED, NWU, SWU, NEU, SEU
 * 12 edges: WD, ED, WU, EU, ND, SD, NU, SU, NE, NW, SE, SW
 * 6 faces: N, S, W, E, D, U
 Input Parameters:
 .  u - input pointer
 
 Output Parameter:
 .  u - output pointer
 .  ctx - optional user-defined context, as set by SNESSetFunction()
 */
PetscErrorCode FormBCu(PetscReal ****u, void*ctx)
{
	PetscErrorCode ierr;
	AppCtx         *user = (AppCtx*)ctx;
	dvi            d = user->d;
	DM             da = user->da;
	PetscInt       mx = user->mx, my = user->my, mz = user->mz;
	PetscReal      n0=user->n0;
	PetscReal      *x=user->x,*y=user->y,*z=user->z;
	PetscReal      L=user->L;
	PetscReal      Zmin=user->outZmin,Z;
	PetscInt       bcType = user->bcType;
	
	PetscInt       i,j,k,l,m,xs,ys,zs,xm,ym,zm;
	PetscInt       rank;
	coeff2         dh;
	coeff3         D1,D2;
	stencil        id;
	
	PetscFunctionBegin;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	// Get local grid boundaries
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	
	
	for (k=zs; k<zs+zm; k++) {
		Z = Zmin + z[k]*L;
		for (j=ys; j<ys+ym; j++) {
			if(xs==0) {		// S-boundary //
				id.x[0] = -1;
				id.x[1] = 0;
				id.x[2] = 1;
				dh.x[0] =    x[id.x[2]] - x[id.x[1]] ;
				dh.x[1] = 2*(x[id.x[2]] - x[id.x[1]]);
				D1.x[0] = - 1.0/dh.x[0] - 1.0/dh.x[1]; 
				D1.x[1] =   1.0/dh.x[0] - 1.0/(dh.x[0]-dh.x[1]); 
				D1.x[2] =   1.0/dh.x[1] + 1.0/(dh.x[0]-dh.x[1]);
				D2.x[0] =   2.0/(dh.x[0] * dh.x[1]); 
				D2.x[1] =   2.0/(dh.x[0] * (dh.x[0]-dh.x[1])); 
				D2.x[2] =   2.0/(dh.x[1] * (dh.x[1]-dh.x[0]));
				
				if (bcType==0) { // 1st order Neumann BC
					for (l=0; l<19; l++) u[k][j][id.x[0]][l] = u[k][j][id.x[1]][l];
				}
				if (bcType==1) { // 2nd Order Neumann BC
					for (l=0; l<19; l++) u[k][j][id.x[0]][l] = -D1.x[1]/D1.x[0]*u[k][j][id.x[1]][l]-D1.x[2]/D1.x[0]*u[k][j][id.x[2]][l];
				}
				if (bcType==2) { // Continuous 2nd derivative BC
					for (l=0; l<19; l++) u[k][j][id.x[0]][l] = -D2.x[1]/D2.x[0]*u[k][j][id.x[1]][l]-D2.x[2]/D2.x[0]*u[k][j][id.x[2]][l];
				}
				if (bcType==10) { // Combination of zero'ed 1st and 2nd derivatives
					for (l=0; l<3; l++) {
						u[k][j][id.x[0]][d.ni[l]] = -D1.x[1]/D1.x[0]*u[k][j][id.x[1]][d.ni[l]]-D1.x[2]/D1.x[0]*u[k][j][id.x[2]][d.ni[l]];
						u[k][j][id.x[0]][d.pi[l]] = -D1.x[1]/D1.x[0]*u[k][j][id.x[1]][d.pi[l]]-D1.x[2]/D1.x[0]*u[k][j][id.x[2]][d.pi[l]];
					}
					u[k][j][id.x[0]][d.pe] = -D1.x[1]/D1.x[0]*u[k][j][id.x[1]][d.pe]-D1.x[2]/D1.x[0]*u[k][j][id.x[2]][d.pe];
					for (m=0; m<3; m++) {
						for (l=0; l<3; l++) {
							u[k][j][id.x[0]][d.vi[l][m]] = -D2.x[1]/D2.x[0]*u[k][j][id.x[1]][d.vi[l][m]]-D2.x[2]/D2.x[0]*u[k][j][id.x[2]][d.vi[l][m]];
						}
						u[k][j][id.x[0]][d.B[m]] = -D1.x[1]/D1.x[0]*u[k][j][id.x[1]][d.B[m]]-D1.x[2]/D1.x[0]*u[k][j][id.x[2]][d.B[m]];
					}
				}
				if (bcType==11) { // Combination of DiRichlet and Neumann
					for (l=0; l< 3; l++) u[k][j][id.x[0]][l] = Interpolate(user->RefPart, l ,Z, lin_flat)*Interpolate(user->RefProf, 7 ,Z, lin_exp)/n0; 
					for (l=3; l<19; l++) u[k][j][id.x[0]][l] = u[k][j][id.x[1]][l];
				}
			}
			if(xs+xm==mx) {		// N-boundary //
				id.x[0] = mx;
				id.x[1] = mx-1;
				id.x[2] = mx-2;
				dh.x[0] =    x[id.x[2]] - x[id.x[1]] ;
				dh.x[1] = 2*(x[id.x[2]] - x[id.x[1]]);
				D1.x[0] = - 1.0/dh.x[0] - 1.0/dh.x[1]; 
				D1.x[1] =   1.0/dh.x[0] - 1.0/(dh.x[0]-dh.x[1]); 
				D1.x[2] =   1.0/dh.x[1] + 1.0/(dh.x[0]-dh.x[1]);
				D2.x[0] =   2.0/(dh.x[0] * dh.x[1]); 
				D2.x[1] =   2.0/(dh.x[0] * (dh.x[0]-dh.x[1])); 
				D2.x[2] =   2.0/(dh.x[1] * (dh.x[1]-dh.x[0]));
				
				if (bcType==0) { // 1st order Neumann BC
					for (l=0; l<19; l++) u[k][j][id.x[0]][l] = u[k][j][id.x[1]][l];
				}
				if (bcType==1) { // 2nd Order Neumann BC
					for (l=0; l<19; l++) u[k][j][id.x[0]][l] = -D1.x[1]/D1.x[0]*u[k][j][id.x[1]][l]-D1.x[2]/D1.x[0]*u[k][j][id.x[2]][l];
				}
				if (bcType==2) { // Continuous 2nd derivative BC
					for (l=0; l<19; l++) u[k][j][id.x[0]][l] = -D2.x[1]/D2.x[0]*u[k][j][id.x[1]][l]-D2.x[2]/D2.x[0]*u[k][j][id.x[2]][l];
				}
				if (bcType==10) { // Combination of zero'ed 1st and 2nd derivatives
					for (l=0; l<3; l++) {
						u[k][j][id.x[0]][d.ni[l]] = -D1.x[1]/D1.x[0]*u[k][j][id.x[1]][d.ni[l]]-D1.x[2]/D1.x[0]*u[k][j][id.x[2]][d.ni[l]];
						u[k][j][id.x[0]][d.pi[l]] = -D1.x[1]/D1.x[0]*u[k][j][id.x[1]][d.pi[l]]-D1.x[2]/D1.x[0]*u[k][j][id.x[2]][d.pi[l]];
					}
					u[k][j][id.x[0]][d.pe] = -D1.x[1]/D1.x[0]*u[k][j][id.x[1]][d.pe]-D1.x[2]/D1.x[0]*u[k][j][id.x[2]][d.pe];
					for (m=0; m<3; m++) {
						for (l=0; l<3; l++) {
							u[k][j][id.x[0]][d.vi[l][m]] = -D2.x[1]/D2.x[0]*u[k][j][id.x[1]][d.vi[l][m]]-D2.x[2]/D2.x[0]*u[k][j][id.x[2]][d.vi[l][m]];
						}
						u[k][j][id.x[0]][d.B[m]] = -D1.x[1]/D1.x[0]*u[k][j][id.x[1]][d.B[m]]-D1.x[2]/D1.x[0]*u[k][j][id.x[2]][d.B[m]];
					}
				}
				if (bcType==11) { // Combination of DiRichlet and Neumann
					for (l=0; l< 3; l++) u[k][j][id.x[0]][l] = Interpolate(user->RefPart, l ,Z, lin_flat)*Interpolate(user->RefProf, 7 ,Z, lin_exp)/n0; 
					for (l=3; l<19; l++) u[k][j][id.x[0]][l] = u[k][j][id.x[1]][l];
				}
			}
		}
	}
	for (i=xs; i<xs+xm; i++) {
		for (k=zs; k<zs+zm; k++) {
			Z = Zmin + z[k]*L;
			if(ys==0) {                     // W-boundary //
				id.y[0] = -1;
				id.y[1] = 0;
				id.y[2] = 1;
				dh.y[0] =    y[id.y[2]] - y[id.y[1]] ;
				dh.y[1] = 2*(y[id.y[2]] - y[id.y[1]]);
				D1.y[0] = - 1.0/dh.y[0] - 1.0/dh.y[1]; 
				D1.y[1] =   1.0/dh.y[0] - 1.0/(dh.y[0]-dh.y[1]); 
				D1.y[2] =   1.0/dh.y[1] + 1.0/(dh.y[0]-dh.y[1]);
				D2.y[0] =   2.0/(dh.y[0] * dh.y[1]); 
				D2.y[1] =   2.0/(dh.y[0] * (dh.y[0]-dh.y[1])); 
				D2.y[2] =   2.0/(dh.y[1] * (dh.y[1]-dh.y[0]));
				
				if (bcType==0) { // 1st order Neumann BC
					for (l=0; l<19; l++) u[k][id.y[0]][i][l] = u[k][id.y[1]][i][l];
				}
				if (bcType==1) { // 2nd Order Neumann BC
					for (l=0; l<19; l++) u[k][id.y[0]][i][l] = -D1.y[1]/D1.y[0]*u[k][id.y[1]][i][l]-D1.y[2]/D1.y[0]*u[k][id.y[2]][i][l];
				}
				if (bcType==2) { // Continuous 2nd derivative BC
					for (l=0; l<19; l++) u[k][id.y[0]][i][l] = -D2.y[1]/D2.y[0]*u[k][id.y[1]][i][l]-D2.y[2]/D2.y[0]*u[k][id.y[2]][i][l];
				}
				if (bcType==10) { // Combination of zero'ed 1st and 2nd derivatives
					for (l=0; l<3; l++) {
						u[k][id.y[0]][i][d.ni[l]] = -D1.y[1]/D1.y[0]*u[k][id.y[1]][i][d.ni[l]]-D1.y[2]/D1.y[0]*u[k][id.y[2]][i][d.ni[l]];
						u[k][id.y[0]][i][d.pi[l]] = -D1.y[1]/D1.y[0]*u[k][id.y[1]][i][d.pi[l]]-D1.y[2]/D1.y[0]*u[k][id.y[2]][i][d.pi[l]];
					}
					u[k][id.y[0]][i][d.pe] = -D1.y[1]/D1.y[0]*u[k][id.y[1]][i][d.pe]-D1.y[2]/D1.y[0]*u[k][id.y[2]][i][d.pe];
					for (m=0; m<3; m++) {
						for (l=0; l<3; l++) {
							u[k][id.y[0]][i][d.vi[l][m]] = -D2.y[1]/D2.y[0]*u[k][id.y[1]][i][d.vi[l][m]]-D2.y[2]/D2.y[0]*u[k][id.y[2]][i][d.vi[l][m]];
						}
						u[k][id.y[0]][i][d.B[m]] = -D1.y[1]/D1.y[0]*u[k][id.y[1]][i][d.B[m]]-D1.y[2]/D1.y[0]*u[k][id.y[2]][i][d.B[m]];
					}
				}
				if (bcType==11) { // Combination of DiRichlet and Neumann
					for (l=0; l< 3; l++) u[k][id.y[0]][i][l] = Interpolate(user->RefPart, l ,Z, lin_flat)*Interpolate(user->RefProf, 7 ,Z, lin_exp)/n0; 
					for (l=3; l<19; l++) u[k][id.y[0]][i][l] = u[k][id.y[1]][i][l] ;
				}
			}
			if(ys+ym==my) {                 // E-boundary //
				id.y[0] = my;
				id.y[1] = my-1;
				id.y[2] = my-2;
				dh.y[0] =    y[id.y[2]] - y[id.y[1]] ;
				dh.y[1] = 2*(y[id.y[2]] - y[id.y[1]]);
				D1.y[0] = - 1.0/dh.y[0] - 1.0/dh.y[1]; 
				D1.y[1] =   1.0/dh.y[0] - 1.0/(dh.y[0]-dh.y[1]); 
				D1.y[2] =   1.0/dh.y[1] + 1.0/(dh.y[0]-dh.y[1]);
				D2.y[0] =   2.0/(dh.y[0] * dh.y[1]); 
				D2.y[1] =   2.0/(dh.y[0] * (dh.y[0]-dh.y[1])); 
				D2.y[2] =   2.0/(dh.y[1] * (dh.y[1]-dh.y[0]));
				
				if (bcType==0) { // 1st order Neumann BC
					for (l=0; l<19; l++) u[k][id.y[0]][i][l] = u[k][id.y[1]][i][l];
				}
				if (bcType==1) { // 2nd Order Neumann BC
					for (l=0; l<19; l++) u[k][id.y[0]][i][l] = -D1.y[1]/D1.y[0]*u[k][id.y[1]][i][l]-D1.y[2]/D1.y[0]*u[k][id.y[2]][i][l];
				}
				if (bcType==2) { // Continuous 2nd derivative BC
					for (l=0; l<19; l++) u[k][id.y[0]][i][l] = -D2.y[1]/D2.y[0]*u[k][id.y[1]][i][l]-D2.y[2]/D2.y[0]*u[k][id.y[2]][i][l];
				}
				if (bcType==10) { // Combination of zero'ed 1st and 2nd derivatives
					for (l=0; l<3; l++) {
						u[k][id.y[0]][i][d.ni[l]] = -D1.y[1]/D1.y[0]*u[k][id.y[1]][i][d.ni[l]]-D1.y[2]/D1.y[0]*u[k][id.y[2]][i][d.ni[l]];
						u[k][id.y[0]][i][d.pi[l]] = -D1.y[1]/D1.y[0]*u[k][id.y[1]][i][d.pi[l]]-D1.y[2]/D1.y[0]*u[k][id.y[2]][i][d.pi[l]];
					}
					u[k][id.y[0]][i][d.pe] = -D1.y[1]/D1.y[0]*u[k][id.y[1]][i][d.pe]-D1.y[2]/D1.y[0]*u[k][id.y[2]][i][d.pe];
					for (m=0; m<3; m++) {
						for (l=0; l<3; l++) {
							u[k][id.y[0]][i][d.vi[l][m]] = -D2.y[1]/D2.y[0]*u[k][id.y[1]][i][d.vi[l][m]]-D2.y[2]/D2.y[0]*u[k][id.y[2]][i][d.vi[l][m]];
						}
						u[k][id.y[0]][i][d.B[m]] = -D1.y[1]/D1.y[0]*u[k][id.y[1]][i][d.B[m]]-D1.y[2]/D1.y[0]*u[k][id.y[2]][i][d.B[m]];
					}
				}
				if (bcType==11) { // Combination of DiRichlet and Neumann
					for (l=0; l< 3; l++) u[k][id.y[0]][i][l] = Interpolate(user->RefPart, l ,Z, lin_flat)*Interpolate(user->RefProf, 7 ,Z, lin_exp)/n0; 
					for (l=3; l<19; l++) u[k][id.y[0]][i][l] = u[k][id.y[1]][i][l] ;
				}
			}
		}
	}
	for (j=ys; j<ys+ym; j++) {
		for (i=xs; i<xs+xm; i++) {
			if(zs==0) {                     // B-boundary //
				Z = Zmin + z[zs]*L;
				
				id.z[0] = -1;
				id.z[1] = 0;
				id.z[2] = 1;
				dh.z[0] =    z[id.z[2]] - z[id.z[1]] ;
				dh.z[1] = 2*(z[id.z[2]] - z[id.z[1]]);
				D1.z[0] = - 1.0/dh.z[0] - 1.0/dh.z[1]; 
				D1.z[1] =   1.0/dh.z[0] - 1.0/(dh.z[0]-dh.z[1]); 
				D1.z[2] =   1.0/dh.z[1] + 1.0/(dh.z[0]-dh.z[1]);
				D2.z[0] =   2.0/(dh.z[0] * dh.z[1]); 
				D2.z[1] =   2.0/(dh.z[0] * (dh.z[0]-dh.z[1])); 
				D2.z[2] =   2.0/(dh.z[1] * (dh.z[1]-dh.z[0]));
				
				if (bcType==0) { // 1st order Neumann BC
					for (l=0; l<19; l++){
						u[id.z[0]][j][i][l] = u[id.z[1]][j][i][l];
					}
				}
				if (bcType==1) { // 2nd Order Neumann BC
					for (l=0; l<19; l++) u[id.z[0]][j][i][l] = -D1.z[1]/D1.z[0]*u[id.z[1]][j][i][l]-D1.z[2]/D1.z[0]*u[id.z[2]][j][i][l];
				}
				if (bcType==2) { // Continuous 2nd derivative BC
					for (l=0; l<19; l++) u[id.z[0]][j][i][l] = -D2.z[1]/D2.z[0]*u[id.z[1]][j][i][l]-D2.z[2]/D2.z[0]*u[id.z[2]][j][i][l];
				}
				if (bcType==10) { // Combination of zero'ed 1st and 2nd derivatives
					for (l=0; l<3; l++) {
						u[id.z[0]][j][i][d.ni[l]] = -D1.z[1]/D1.z[0]*u[id.z[1]][j][i][d.ni[l]]-D1.z[2]/D1.z[0]*u[id.z[2]][j][i][d.ni[l]];
						u[id.z[0]][j][i][d.pi[l]] = -D1.z[1]/D1.z[0]*u[id.z[1]][j][i][d.pi[l]]-D1.z[2]/D1.z[0]*u[id.z[2]][j][i][d.pi[l]];
					}
					u[id.z[0]][j][i][d.pe] = -D1.z[1]/D1.z[0]*u[id.z[1]][j][i][d.pe]-D1.z[2]/D1.z[0]*u[id.z[2]][j][i][d.pe];
					for (m=0; m<3; m++) {
						for (l=0; l<3; l++) {
							u[id.z[0]][j][i][d.vi[l][m]] = -D2.z[1]/D2.z[0]*u[id.z[1]][j][i][d.vi[l][m]]-D2.z[2]/D2.z[0]*u[id.z[2]][j][i][d.vi[l][m]];
						}
						u[id.z[0]][j][i][d.B[m]] = -D1.z[1]/D1.z[0]*u[id.z[1]][j][i][d.B[m]]-D1.z[2]/D1.z[0]*u[id.z[2]][j][i][d.B[m]];
					}
				}
				if (bcType==11) { // Combination of DiRichlet and Neumann
					for (l=0; l< 3; l++) u[id.z[0]][j][i][l] = Interpolate(user->RefPart, l ,Z, lin_flat)*Interpolate(user->RefProf, 7 ,Z, lin_exp)/n0; 
					for (l=3; l<19; l++) u[id.z[0]][j][i][l] = u[id.z[1]][j][i][l] ;
				}
			}
			if(zs+zm==mz) {                 // T-boundary //
				Z = Zmin + z[mz-1]*L;
				
				id.z[0] = mz;
				id.z[1] = mz-1;
				id.z[2] = mz-2;
				dh.z[0] =    z[id.z[2]] - z[id.z[1]] ;
				dh.z[1] = 2*(z[id.z[2]] - z[id.z[1]]);
				D1.z[0] = - 1.0/dh.z[0] - 1.0/dh.z[1]; 
				D1.z[1] =   1.0/dh.z[0] - 1.0/(dh.z[0]-dh.z[1]); 
				D1.z[2] =   1.0/dh.z[1] + 1.0/(dh.z[0]-dh.z[1]);
				D2.z[0] =   2.0/(dh.z[0] * dh.z[1]); 
				D2.z[1] =   2.0/(dh.z[0] * (dh.z[0]-dh.z[1])); 
				D2.z[2] =   2.0/(dh.z[1] * (dh.z[1]-dh.z[0]));
				
				if (bcType==0) { // 1st order Neumann BC
					for (l=0; l<19; l++) u[id.z[0]][j][i][l] = u[id.z[1]][j][i][l];
				}
				if (bcType==1) { // 2nd Order Neumann BC
					for (l=0; l<19; l++) u[id.z[0]][j][i][l] = -D1.z[1]/D1.z[0]*u[id.z[1]][j][i][l]-D1.z[2]/D1.z[0]*u[id.z[2]][j][i][l];
				}
				if (bcType==2) { // Continuous 2nd derivative BC
					for (l=0; l<19; l++) u[id.z[0]][j][i][l] = -D2.z[1]/D2.z[0]*u[id.z[1]][j][i][l]-D2.z[2]/D2.z[0]*u[id.z[2]][j][i][l];
				}
				if (bcType==10) { // Combination of zero'ed 1st and 2nd derivatives
					for (l=0; l<3; l++) {
						u[id.z[0]][j][i][d.ni[l]] = -D1.z[1]/D1.z[0]*u[id.z[1]][j][i][d.ni[l]]-D1.z[2]/D1.z[0]*u[id.z[2]][j][i][d.ni[l]];
						u[id.z[0]][j][i][d.pi[l]] = -D1.z[1]/D1.z[0]*u[id.z[1]][j][i][d.pi[l]]-D1.z[2]/D1.z[0]*u[id.z[2]][j][i][d.pi[l]];
					}
					u[id.z[0]][j][i][d.pe] = -D1.z[1]/D1.z[0]*u[id.z[1]][j][i][d.pe]-D1.z[2]/D1.z[0]*u[id.z[2]][j][i][d.pe];
					for (m=0; m<3; m++) {
						for (l=0; l<3; l++) {
							u[id.z[0]][j][i][d.vi[l][m]] = -D2.z[1]/D2.z[0]*u[id.z[1]][j][i][d.vi[l][m]]-D2.z[2]/D2.z[0]*u[id.z[2]][j][i][d.vi[l][m]];
						}
						u[id.z[0]][j][i][d.B[m]] = -D1.z[1]/D1.z[0]*u[id.z[1]][j][i][d.B[m]]-D1.z[2]/D1.z[0]*u[id.z[2]][j][i][d.B[m]];
					}
				}
				if (bcType==11) { // Combination of DiRichlet and Neumann
					for (l=0; l< 3; l++) u[id.z[0]][j][i][l] = Interpolate(user->RefPart, l ,Z, lin_flat)*Interpolate(user->RefProf, 7 ,Z, lin_exp)/n0; 
					for (l=3; l<19; l++) u[id.z[0]][j][i][l] = u[id.z[1]][j][i][l] ;
				}
			}
		}
	}
	PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormBCv"
/* 
 FormBCv - Evaluates intermediate function @ BC (outside the domain), F(u).
 * T: Top, B: Bottom, N: North, S: South, E: East, W: West
 * 8 corners: NWD, SWD, NED, SED, NWU, SWU, NEU, SEU
 * 12 edges: WD, ED, WU, EU, ND, SD, NU, SU, NE, NW, SE, SW
 * 6 faces: N, S, W, E, D, U
 Input Parameters:
 .  v - input pointer
 
 Output Parameter:
 .  v - output pointer
 .  ctx - optional user-defined context, as set by SNESSetFunction()
 */

PetscErrorCode FormBCv(PetscReal ****v, void *ctx)
{
	PetscErrorCode ierr;
	AppCtx         *user = (AppCtx*)ctx;
	DM             da = (DM)user->da;
	svi            s = user->s;
	
	PetscInt       mx = user->mx, my = user->my, mz = user->mz;
	PetscInt       rank;
	PetscInt       i,j,k,l,m,xs,ys,zs,xm,ym,zm;
	
	PetscReal      *x=user->x,*y=user->y,*z=user->z;
	stencil        id;
	PetscInt       bcType = user->bcType;
	coeff3         D1,D2;
	coeff2         dh;
	
	PetscFunctionBegin;
	
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	id.x[0] = 0; id.y[0] = 0; id.z[0] = 0;
	id.x[1] = 0; id.y[1] = 0; id.z[1] = 0;
	id.x[2] = 0; id.y[2] = 0; id.z[2] = 0;
	
	dh.x[0] = 0; dh.y[0] = 0; dh.z[0] = 0;
	dh.x[1] = 0; dh.y[1] = 0; dh.z[1] = 0;
	
	// Get local grid boundaries
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	
	for (j=ys; j<ys+ym; j++) {
		for (i=xs; i<xs+xm; i++) {
			if(zs==0) {                     // B-boundary //
				id.z[0] = -1;
				id.z[1] = 0;
				id.z[2] = 1;
				dh.z[0] =    z[id.z[2]] - z[id.z[1]] ;
				dh.z[1] = 2*(z[id.z[2]] - z[id.z[1]]);
				D1.z[0] = - 1.0/dh.z[0] - 1.0/dh.z[1]; 
				D1.z[1] =   1.0/dh.z[0] - 1.0/(dh.z[0]-dh.z[1]); 
				D1.z[2] =   1.0/dh.z[1] + 1.0/(dh.z[0]-dh.z[1]);
				D2.z[0] =   2.0/(dh.z[0] * dh.z[1]); 
				D2.z[1] =   2.0/(dh.z[0] * (dh.z[0]-dh.z[1])); 
				D2.z[2] =   2.0/(dh.z[1] * (dh.z[1]-dh.z[0]));
				
				if (bcType==0 || bcType==11) { // 1st order Neumann BC
					for (l=0; l<9; l++) v[id.z[0]][j][i][l] = v[id.z[1]][j][i][l];
				}
				if (bcType==1) { // 2nd Order Neumann BC
					for (l=0; l<9; l++) v[id.z[0]][j][i][l] = -D1.z[1]/D1.z[0]*v[id.z[1]][j][i][l]-D1.z[2]/D1.z[0]*v[id.z[2]][j][i][l];
				}
				if (bcType==2) { // Continuous 2nd derivative BC
					for (l=0; l<9; l++) v[id.z[0]][j][i][l] = -D2.z[1]/D2.z[0]*v[id.z[1]][j][i][l]-D2.z[2]/D2.z[0]*v[id.z[2]][j][i][l];
				}
				if (bcType==10) { // Combination of zero'ed 1st and 2nd derivatives
					for (m=0; m<3; m++) {
						v[id.z[0]][j][i][s.ve[m]] = -D2.z[1]/D2.z[0]*v[id.z[1]][j][i][s.ve[m]]-D2.z[2]/D2.z[0]*v[id.z[2]][j][i][s.ve[m]];
						v[id.z[0]][j][i][s.E[m]]  = -D1.z[1]/D1.z[0]*v[id.z[1]][j][i][s.E[m]] -D1.z[2]/D1.z[0]*v[id.z[2]][j][i][s.E[m]];
						v[id.z[0]][j][i][s.J[m]]  = -D1.z[1]/D1.z[0]*v[id.z[1]][j][i][s.J[m]] -D1.z[2]/D1.z[0]*v[id.z[2]][j][i][s.J[m]];
					}
				}
			}
			if(zs+zm==mz) {                 // T-boundary //
				id.z[0] = mz;
				id.z[1] = mz-1;
				id.z[2] = mz-2;
				dh.z[0] =    z[id.z[2]] - z[id.z[1]] ;
				dh.z[1] = 2*(z[id.z[2]] - z[id.z[1]]);
				D1.z[0] = - 1.0/dh.z[0] - 1.0/dh.z[1]; 
				D1.z[1] =   1.0/dh.z[0] - 1.0/(dh.z[0]-dh.z[1]); 
				D1.z[2] =   1.0/dh.z[1] + 1.0/(dh.z[0]-dh.z[1]);
				D2.z[0] =   2.0/(dh.z[0] * dh.z[1]); 
				D2.z[1] =   2.0/(dh.z[0] * (dh.z[0]-dh.z[1])); 
				D2.z[2] =   2.0/(dh.z[1] * (dh.z[1]-dh.z[0]));
				
				if (bcType==0 || bcType==11) { // 1st order Neumann BC
					for (l=0; l<9; l++) v[id.z[0]][j][i][l] = v[id.z[1]][j][i][l];
				}
				if (bcType==1) { // 2nd Order Neumann BC
					for (l=0; l<9; l++) v[id.z[0]][j][i][l] = -D1.z[1]/D1.z[0]*v[id.z[1]][j][i][l]-D1.z[2]/D1.z[0]*v[id.z[2]][j][i][l];
				}
				if (bcType==2) { // Continuous 2nd derivative BC
					for (l=0; l<9; l++) v[id.z[0]][j][i][l] = -D2.z[1]/D2.z[0]*v[id.z[1]][j][i][l]-D2.z[2]/D2.z[0]*v[id.z[2]][j][i][l];
				}
				if (bcType==10) { // Combination of zero'ed 1st and 2nd derivatives
					for (m=0; m<3; m++) {
						v[id.z[0]][j][i][s.ve[m]] = -D2.z[1]/D2.z[0]*v[id.z[1]][j][i][s.ve[m]]-D2.z[2]/D2.z[0]*v[id.z[2]][j][i][s.ve[m]];
						v[id.z[0]][j][i][s.E[m]]  = -D1.z[1]/D1.z[0]*v[id.z[1]][j][i][s.E[m]] -D1.z[2]/D1.z[0]*v[id.z[2]][j][i][s.E[m]];
						v[id.z[0]][j][i][s.J[m]]  = -D1.z[1]/D1.z[0]*v[id.z[1]][j][i][s.J[m]] -D1.z[2]/D1.z[0]*v[id.z[2]][j][i][s.J[m]];
					}
				}
			}
		}
	}
	
	for (k=zs; k<zs+zm; k++) {
		for (j=ys; j<ys+ym; j++) {
			if(xs==0) {                     // W-boundary //
				id.x[0] = -1;
				id.x[1] = 0;
				id.x[2] = 1;
				dh.x[0] =    x[id.x[2]] - x[id.x[1]] ;
				dh.x[1] = 2*(x[id.x[2]] - x[id.x[1]]);
				D1.x[0] = - 1.0/dh.x[0] - 1.0/dh.x[1]; 
				D1.x[1] =   1.0/dh.x[0] - 1.0/(dh.x[0]-dh.x[1]); 
				D1.x[2] =   1.0/dh.x[1] + 1.0/(dh.x[0]-dh.x[1]);
				D2.x[0] =   2.0/(dh.x[0] * dh.x[1]); 
				D2.x[1] =   2.0/(dh.x[0] * (dh.x[0]-dh.x[1])); 
				D2.x[2] =   2.0/(dh.x[1] * (dh.x[1]-dh.x[0]));
				
				if (bcType==0 || bcType==11) { // 1st order Neumann BC
					for (l=0; l<9; l++) v[k][j][id.x[0]][l] = v[k][j][id.x[1]][l];
				}
				if (bcType==1) { // 2nd Order Neumann BC
					for (l=0; l<9; l++) v[k][j][id.x[0]][l] = -D1.x[1]/D1.x[0]*v[k][j][id.x[1]][l]-D1.x[2]/D1.x[0]*v[k][j][id.x[2]][l];
				}
				if (bcType==2) { // Continuous 2nd derivative BC
					for (l=0; l<9; l++) v[k][j][id.x[0]][l] = -D2.x[1]/D2.x[0]*v[k][j][id.x[1]][l]-D2.x[2]/D2.x[0]*v[k][j][id.x[2]][l];
				}
				if (bcType==10) { // Combination of zero'ed 1st and 2nd derivatives
					for (m=0; m<3; m++) {
						v[k][j][id.x[0]][s.ve[m]] = -D2.x[1]/D2.x[0]*v[k][j][id.x[1]][s.ve[m]]-D2.x[2]/D2.x[0]*v[k][j][id.x[2]][s.ve[m]];
						v[k][j][id.x[0]][s.E[m]]  = -D1.x[1]/D1.x[0]*v[k][j][id.x[1]][s.E[m]] -D1.x[2]/D1.x[0]*v[k][j][id.x[2]][s.E[m]];
						v[k][j][id.x[0]][s.J[m]]  = -D1.x[1]/D1.x[0]*v[k][j][id.x[1]][s.J[m]] -D1.x[2]/D1.x[0]*v[k][j][id.x[2]][s.J[m]];
					}
				}
			}
			if(xs+xm==mx) {                 // E-boundary //
				id.x[0] = mx;
				id.x[1] = mx-1;
				id.x[2] = mx-2;
				dh.x[0] =    x[id.x[2]] - x[id.x[1]] ;
				dh.x[1] = 2*(x[id.x[2]] - x[id.x[1]]);
				D1.x[0] = - 1.0/dh.x[0] - 1.0/dh.x[1]; 
				D1.x[1] =   1.0/dh.x[0] - 1.0/(dh.x[0]-dh.x[1]); 
				D1.x[2] =   1.0/dh.x[1] + 1.0/(dh.x[0]-dh.x[1]);
				D2.x[0] =   2.0/(dh.x[0] * dh.x[1]); 
				D2.x[1] =   2.0/(dh.x[0] * (dh.x[0]-dh.x[1])); 
				D2.x[2] =   2.0/(dh.x[1] * (dh.x[1]-dh.x[0]));
				
				if (bcType==0 || bcType==11) { // 1st order Neumann BC
					for (l=0; l<9; l++) v[k][j][id.x[0]][l] = v[k][j][id.x[1]][l];
				}
				if (bcType==1) { // 2nd Order Neumann BC
					for (l=0; l<9; l++) v[k][j][id.x[0]][l] = -D1.x[1]/D1.x[0]*v[k][j][id.x[1]][l]-D1.x[2]/D1.x[0]*v[k][j][id.x[2]][l];
				}
				if (bcType==2) { // Continuous 2nd derivative BC
					for (l=0; l<9; l++) v[k][j][id.x[0]][l] = -D2.x[1]/D2.x[0]*v[k][j][id.x[1]][l]-D2.x[2]/D2.x[0]*v[k][j][id.x[2]][l];
				}
				if (bcType==10) { // Combination of zero'ed 1st and 2nd derivatives
					for (m=0; m<3; m++) {
						v[k][j][id.x[0]][s.ve[m]] = -D2.x[1]/D2.x[0]*v[k][j][id.x[1]][s.ve[m]]-D2.x[2]/D2.x[0]*v[k][j][id.x[2]][s.ve[m]];
						v[k][j][id.x[0]][s.E[m]]  = -D1.x[1]/D1.x[0]*v[k][j][id.x[1]][s.E[m]] -D1.x[2]/D1.x[0]*v[k][j][id.x[2]][s.E[m]];
						v[k][j][id.x[0]][s.J[m]]  = -D1.x[1]/D1.x[0]*v[k][j][id.x[1]][s.J[m]] -D1.x[2]/D1.x[0]*v[k][j][id.x[2]][s.J[m]];
					}
				}
			}
		}
	}
	for (i=xs; i<xs+xm; i++) {
		for (k=zs; k<zs+zm; k++) {
			if(ys==0) {                     // S-boundary //
				id.y[0] = -1;
				id.y[1] = 0;
				id.y[2] = 1;
				dh.y[0] =    y[id.y[2]] - y[id.y[1]] ;
				dh.y[1] = 2*(y[id.y[2]] - y[id.y[1]]);
				D1.y[0] = - 1.0/dh.y[0] - 1.0/dh.y[1]; 
				D1.y[1] =   1.0/dh.y[0] - 1.0/(dh.y[0]-dh.y[1]); 
				D1.y[2] =   1.0/dh.y[1] + 1.0/(dh.y[0]-dh.y[1]);
				D2.y[0] =   2.0/(dh.y[0] * dh.y[1]); 
				D2.y[1] =   2.0/(dh.y[0] * (dh.y[0]-dh.y[1])); 
				D2.y[2] =   2.0/(dh.y[1] * (dh.y[1]-dh.y[0]));
				
				if (bcType==0 || bcType==11) { // 1st order Neumann BC
					for (l=0; l<9; l++) v[k][id.y[0]][i][l] = v[k][id.y[1]][i][l];
				}
				if (bcType==1) { // 2nd Order Neumann BC
					for (l=0; l<9; l++) v[k][id.y[0]][i][l] = -D1.y[1]/D1.y[0]*v[k][id.y[1]][i][l]-D1.y[2]/D1.y[0]*v[k][id.y[2]][i][l];
				}
				if (bcType==2) { // Continuous 2nd derivative BC
					for (l=0; l<9; l++) v[k][id.y[0]][i][l] = -D2.y[1]/D2.y[0]*v[k][id.y[1]][i][l]-D2.y[2]/D2.y[0]*v[k][id.y[2]][i][l];
				}
				if (bcType==10) { // Combination of zero'ed 1st and 2nd derivatives
					for (m=0; m<3; m++) {
						v[k][id.y[0]][i][s.ve[m]] = -D2.y[1]/D2.y[0]*v[k][id.y[1]][i][s.ve[m]]-D2.y[2]/D2.y[0]*v[k][id.y[2]][i][s.ve[m]];
						v[k][id.y[0]][i][s.E[m]]  = -D1.y[1]/D1.y[0]*v[k][id.y[1]][i][s.E[m]] -D1.y[2]/D1.y[0]*v[k][id.y[2]][i][s.E[m]];
						v[k][id.y[0]][i][s.J[m]]  = -D1.y[1]/D1.y[0]*v[k][id.y[1]][i][s.J[m]] -D1.y[2]/D1.y[0]*v[k][id.y[2]][i][s.J[m]];
					}
				}
			}
			if(ys+ym==my) {                 // N-boundary //
				id.y[0] = my;
				id.y[1] = my-1;
				id.y[2] = my-2;
				dh.y[0] =    y[id.y[2]] - y[id.y[1]] ;
				dh.y[1] = 2*(y[id.y[2]] - y[id.y[1]]);
				D1.y[0] = - 1.0/dh.y[0] - 1.0/dh.y[1]; 
				D1.y[1] =   1.0/dh.y[0] - 1.0/(dh.y[0]-dh.y[1]); 
				D1.y[2] =   1.0/dh.y[1] + 1.0/(dh.y[0]-dh.y[1]);
				D2.y[0] =   2.0/(dh.y[0] * dh.y[1]); 
				D2.y[1] =   2.0/(dh.y[0] * (dh.y[0]-dh.y[1])); 
				D2.y[2] =   2.0/(dh.y[1] * (dh.y[1]-dh.y[0]));
				
				if (bcType==0 || bcType==11) { // 1st order Neumann BC
					for (l=0; l<9; l++) v[k][id.y[0]][i][l] = v[k][id.y[1]][i][l];
				}
				if (bcType==1) { // 2nd Order Neumann BC
					for (l=0; l<9; l++) v[k][id.y[0]][i][l] = -D1.y[1]/D1.y[0]*v[k][id.y[1]][i][l]-D1.y[2]/D1.y[0]*v[k][id.y[2]][i][l];
				}
				if (bcType==2) { // Continuous 2nd derivative BC
					for (l=0; l<9; l++) v[k][id.y[0]][i][l] = -D2.y[1]/D2.y[0]*v[k][id.y[1]][i][l]-D2.y[2]/D2.y[0]*v[k][id.y[2]][i][l];
				}
				if (bcType==10) { // Combination of zero'ed 1st and 2nd derivatives
					for (m=0; m<3; m++) {
						v[k][id.y[0]][i][s.ve[m]] = -D2.y[1]/D2.y[0]*v[k][id.y[1]][i][s.ve[m]]-D2.y[2]/D2.y[0]*v[k][id.y[2]][i][s.ve[m]];
						v[k][id.y[0]][i][s.E[m]]  = -D1.y[1]/D1.y[0]*v[k][id.y[1]][i][s.E[m]] -D1.y[2]/D1.y[0]*v[k][id.y[2]][i][s.E[m]];
						v[k][id.y[0]][i][s.J[m]]  = -D1.y[1]/D1.y[0]*v[k][id.y[1]][i][s.J[m]] -D1.y[2]/D1.y[0]*v[k][id.y[2]][i][s.J[m]];
					}
				}
			}
		}
	}   
	PetscFunctionReturn(0); 
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "CheckValues"
/* 
 CheckValues - Check the domain for unphysical values of density, pressure, etc..
 
 Input Parameters:
 .  U - input vector
 .  V - input vector
 .  ctx - optional user-defined context, as set by SNESSetFunction()
 
 Output Parameter:
 .  U - output vector
 .  V - output vector
 */
PetscErrorCode CheckValues(PetscReal ****u,PetscReal ****v,void *ctx)
{
	PetscErrorCode ierr;
	AppCtx         *user = (AppCtx*)ctx;
	DM             da = (DM)user->da;
	dvi            d = user->d;
	svi            s = user->s;
	tolerance      eps = user->eps;
	PetscBool      limiters = user->limiters;
	PetscBool      blockers = user->blockers;
	PetscReal      n0 = user->n0, v0 = user->v0, p0 = user->p0;
	PetscInt       rank;
	PetscReal      *z=user->z;
	PetscReal      Z;
	PetscReal      Zmin = user->outZmin;
	PetscReal      L=user->L;
	PetscInt       i,j,k,l,m,xs,ys,zs,xm,ym,zm;
	
	PetscFunctionBegin;
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	
	// Get local grid boundaries
	ierr = DMDAGetGhostCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	
	// Check values
	for (k=zs; k<zs+zm; k++) {
		Z = Zmin + z[k]*L; 
		for (j=ys; j<ys+ym; j++) {
			for (i=xs; i<xs+xm; i++) {
				if (limiters) {
					if(!(u[k][j][i][d.pe]>=eps.p)) {
						u[k][j][i][d.pe]=eps.p;
					}
					assert(u[k][j][i][d.pe]>=eps.p); // check pe>=0
					for (l=0; l<3; l++) {
						if(!(u[k][j][i][d.ni[l]]>=eps.n)) {
							u[k][j][i][d.ni[l]]=eps.n;
						}
						assert(u[k][j][i][d.ni[l]]>=eps.n); // check ni>0
						if(!(u[k][j][i][d.pi[l]]>=eps.p)) {
							u[k][j][i][d.pi[l]]=eps.p;
							assert(u[k][j][i][d.pi[l]]>=eps.p); // check pi>=0
						}
					}
				}
				if (blockers) {
					if(!(u[k][j][i][d.pe]>=0)) PetscPrintf(PETSC_COMM_SELF, "pe=%2.3e @ [%d,%d,%d] (Z=%f km)\n",u[k][j][i][d.pe]*p0, i,j,k, Z*1e-3);
					assert(u[k][j][i][d.pe]>=0); // check pe>=0
					for (l=0; l<3; l++) {
						if(!(u[k][j][i][d.ni[l]]>0)) PetscPrintf(PETSC_COMM_SELF, "ni[%d]=%2.3e @ [%d,%d,%d] (Z=%f km)\n",l,u[k][j][i][d.ni[l]]*n0, i,j,k, Z*1e-3);
						assert(u[k][j][i][d.ni[l]]>0); // check ni>0
						if(!(u[k][j][i][d.pi[l]]>=0)) PetscPrintf(PETSC_COMM_SELF, "pi[%d]=%2.3e @ [%d,%d,%d] (Z=%f km)\n",l,u[k][j][i][d.pi[l]]*p0, i,j,k, Z*1e-3);
						assert(u[k][j][i][d.pi[l]]>=0); // check pi>=0
						for (m=0; m<3; m++) {
							if(!(u[k][j][i][d.vi[l][m]]<=eps.v)) PetscPrintf(PETSC_COMM_SELF, "vi[%d][%d]=%2.3e @ [%d,%d,%d] (Z=%f km)\n",l,m,u[k][j][i][d.vi[l][m]]*v0, i,j,k, Z*1e-3);
							assert(u[k][j][i][d.vi[l][m]]<=eps.v); // check vi<=eps.v
						}
					}
					for (m=0; m<3; m++) {
						if(!(v[k][j][i][s.ve[m]]<=eps.v)) PetscPrintf(PETSC_COMM_SELF, "ve[%d]=%2.3e @ [%d,%d,%d] (Z=%f km)\n",m,v[k][j][i][s.ve[m]]*v0, i,j,k, Z*1e-3);
						assert(v[k][j][i][s.ve[m]]<=eps.v); // check vi<=eps.v
					}
				}
			}
		}
	}
	PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "CalculateFluxes"
/* 
 CalculateFlows - Evaluates Fluxes at BC of the simulation domain.
 
 Input Parameters:
 .  t - current time
 .  u - input pointer
 .  v - input pointer
 
 Output Parameter:
 .  ctx - optional user-defined context, as set by SNESSetFunction()
 */
PetscErrorCode CalculateFluxes(PetscReal t, PetscInt step, PetscReal ****u, PetscReal ****v, void *ctx)
{
	
	PetscErrorCode ierr;
	AppCtx         *user = (AppCtx*)ctx;
	PetscReal      Zmin = user->outZmin;
	PetscReal      n0 = user->n0, v0 = user->v0;
	PetscReal      *x = user->x, *y = user->y, *z = user->z;
	PetscReal      DX,DY,DZ;
	PetscInt       mx = user->mx, my = user->my, mz = user->mz;
	DM             da = user->da;
	PetscReal      L = user->L;
	dvi            d = user->d;
	svi            s = user->s;
	
	FILE           *nFile;
	char           fName[PETSC_MAX_PATH_LEN];
	PetscInt       rank;
	PetscInt       i,j,k,l,m,xs,ys,zs,xm,ym,zm;
	PetscInt       imin,imax,jmin,jmax,kmin,kmax;
	PetscReal      ne=0,ni[3]={0.0,0.0,0.0},ve[3],vi[3][3];
	PetscReal      Fe[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	PetscReal      Fi[3][6] = { {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
	PetscReal      Z;
	
	PetscFunctionBegin;
	
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	
	// Get local grid boundaries
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	
	imin = 0;
	imax = mx-1;
	jmin = 0;
	jmax = my-1;
	kmin = 0;
	kmax = mz-1;
	
	// Calculate the flows through the boundaries of the PLOTTED domain
	PetscPrintf(PETSC_COMM_WORLD,"START FLUX CALCULATIONS.\n");
	// W-E flows
	for (j=ys; j<ys+ym; j++) {
		// Define y-space increment
		if (j==0)         DY = (y[   1]-y[   0])*L    ;
		else if (j==my-1) DY = (y[my-1]-y[my-2])*L    ;
		else              DY = (y[ j+1]-y[ j-1])*L/2.0;
		
		for (k=zs; k<zs+zm; k++) {
			// Define z-space increment
			if (k==0)         DZ = (z[   1]-z[   0])*L    ;
			else if (k==mz-1) DZ = (z[mz-1]-z[mz-2])*L    ;
			else              DZ = (z[ k+1]-z[ k-1])*L/2.0;
			
			// W-flow
			i = imin;
			ne = 0.0;
			for (l=0; l<3; l++) {
				ni[l]    = u[k][j][i][d.ni[l]];
				vi[l][0] = u[k][j][i][d.vi[l][0]];
				ne += ni[l];
				Fi[l][0] -= ni[l]*n0 * vi[l][0]*v0 * DY*DZ;
			}
			ve[0] = v[k][j][i][s.ve[0]];
			Fe[0] -= ne*n0 * ve[0]*v0 * DY*DZ;
			
			// E-flow
			i = imax;
			ne = 0.0;
			for (l=0; l<3; l++) {
				ni[l] = u[k][j][i][d.ni[l]];
				vi[l][0] = u[k][j][i][d.vi[l][0]];
				ne += ni[l];
				Fi[l][1] += ni[l]*n0 * vi[l][0]*v0 * DY*DZ;
			}
			ve[0] = v[k][j][i][s.ve[0]];
			Fe[1] += ne*n0 * ve[0]*v0 * DY*DZ;
		}
	}
	
	// S-N flows
	for (k=zs; k<zs+zm; k++) {
		// Define z-space increment
		if (k==0)         DZ = (z[   1]-z[   0])*L    ;
		else if (k==mz-1) DZ = (z[mz-1]-z[mz-2])*L    ;
		else              DZ = (z[ k+1]-z[ k-1])*L/2.0;
		
		for (i=xs; i<xs+xm; i++) {
			
			// Define x-space increment
			if (i==0)         DX = (x[   1]-x[   0])*L    ;
			else if (i==mx-1) DX = (x[mx-1]-x[mx-2])*L    ;
			else              DX = (x[ i+1]-x[ i-1])*L/2.0;
			
			// S-flow
			j = jmin;
			ne = 0.0;
			for (l=0; l<3; l++) {
				ni[l] = u[k][j][i][d.ni[l]];
				vi[l][1] = u[k][j][i][d.vi[l][1]];
				ne += ni[l];
				Fi[l][2] -= ni[l]*n0 * vi[l][1]*v0 * DZ*DX;
			}
			ve[1] = v[k][j][i][s.ve[1]];
			Fe[2] -= ne*n0 * ve[1]*v0 * DZ*DX;
			
			// N-flow
			j = jmax;
			ne = 0.0;
			for (l=0; l<3; l++) {
				ni[l] = u[k][j][i][d.ni[l]];
				vi[l][1] = u[k][j][i][d.vi[l][1]];
				ne += ni[l];
				Fi[l][3] += ni[l]*n0 * vi[l][1]*v0 * DZ*DX;
			}
			ve[1] = v[k][j][i][s.ve[1]];
			Fe[3] += ne*n0 * ve[1]*v0 * DZ*DX;
			
		}
	}
	
	// B-T flows
	for (i=xs; i<xs+xm; i++) {
		// Define x-space increment
		if (i==0)         DX = (x[   1]-x[   0])*L    ;
		else if (i==mx-1) DX = (x[mx-1]-x[mx-2])*L    ;
		else              DX = (x[ i+1]-x[ i-1])*L/2.0;
		
		for (j=ys; j<ys+ym; j++) {
			// Define y-space increment
			if (j==0)         DY = (y[   1]-y[   0])*L    ;
			else if (j==my-1) DY = (y[my-1]-y[my-2])*L    ;
			else              DY = (y[ j+1]-y[ j-1])*L/2.0;
			
			// B-flow
			k = kmin;
			ne = 0.0;
			for (l=0; l<3; l++) {
				ni[l] = u[k][j][i][d.ni[l]];
				vi[l][2] = u[k][j][i][d.vi[l][2]];
				ne += ni[l]; 
				Fi[l][4] -= ni[l]*n0 * vi[l][2]*v0 * DX*DY;
			}
			ve[2] = v[k][j][i][s.ve[2]];
			Fe[4] -= ne*n0 *ve[2]*v0 * DX*DY;
			
			// T-flow
			k = kmax;
			ne = 0.0;
			for (l=0; l<3; l++) {
				ni[l] = u[k][j][i][d.ni[l]];
				vi[l][2] = u[k][j][i][d.vi[l][2]];
				ne += ni[l];
				Fi[l][5] += ni[l]*n0 * vi[l][2]*v0 * DX*DY;
			}
			ve[2] = v[k][j][i][s.ve[2]];
			Fe[5] += ne*n0 * ve[2]*v0 * DX*DY;
		}
	}
	PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormIntermediateFunction"
/* 
 FormIntermediateFunction - Evaluates nonlinear function, V=G(U).
 
 Input Parameters:
 .  U - input vector
 
 Output Parameter:
 .  V - output vector
 .  ctx - optional user-defined context, as set by SNESSetFunction()
 */
PetscErrorCode FormIntermediateFunction(PetscReal ****u, Vec V, void *ctx)
{
	PetscErrorCode ierr;
	AppCtx         *user = (AppCtx*)ctx;
	DM             da = (DM)user->da;
	DM             db = (DM)user->db;
	svi            s = user->s;
	
	PetscInt       mx = user->mx, my = user->my, mz = user->mz;
	PetscReal      n0 = user->n0, v0 = user->v0, p0 = user->p0, T0 = user->T0;
	PetscReal      un[3] = {user->un[0]/v0, user->un[1]/v0, user->un[2]/v0};       // neutral wind
	PetscInt       rank;
	PetscInt       b,i,j,k,m,xs,ys,zs,xm,ym,zm;
	PetscReal      vi[3][3],ve[3],vb[6][3]; // velocity of i,e
	PetscReal      E[3];  // E-field
	PetscReal      J[3];  // total conduction current
	PetscReal      el_coll[3] = {0.0,0.0,0.0}; // Elastic collision term in the GIVEN direction m  at the CURRENT point (k,j,i)
	PetscReal      nu[6] = {0.0,0.0,0.0,0.0,0.0,0.0}; // electron-neutral collisions
	dvi            d = user->d;
	PetscReal      L = user->L;
	PetscReal      Lz;
	
	PetscBool      vDamping = user->vDamping;
	PetscReal      lambda=user->lambda;
	PetscReal      XL=user->inXmin,XU=user->inXmax;
	PetscReal      YL=user->inYmin,YU=user->inYmax;
	PetscReal      ZL=user->ZL    ,ZU=user->ZU    ;
	PetscReal      Xmin=user->outXmin,Ymin=user->outYmin, Zmin=user->outZmin;
	PetscReal      X,Y,Z;
	
	PetscReal      *x=user->x,*y=user->y,*z=user->z;
	PetscReal      tau = user->tau; 
	PetscReal      ni[3],ne,nn[2]; // density of i,e
	PetscReal      pe;
	PetscReal      Te;
	PetscReal      B[3];
	diff           dB[3],dpe;
	stencil        id;
	PetscReal      ****v;
	PetscInt       bcType = user->bcType;
	coeff3         D1,D2;
	coeff2         dh;
    PetscReal      chem_nuS[7][5];
    PetscReal      E_chem_S[3];

	
	PetscFunctionBegin;
	
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	Lz = user->outZmax - user->outZmin;
	
	ierr = DMDAVecGetArrayDOF(db,V,&v);CHKERRQ(ierr);
	id.x[0] = 0; id.y[0] = 0; id.z[0] = 0;
	id.x[1] = 0; id.y[1] = 0; id.z[1] = 0;
	id.x[2] = 0; id.y[2] = 0; id.z[2] = 0;
	
	dh.x[0] = 0; dh.y[0] = 0; dh.z[0] = 0;
	dh.x[1] = 0; dh.y[1] = 0; dh.z[1] = 0;
	
	// Get local grid boundaries
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
	
	for (k=zs; k<zs+zm; k++) {
		Z = Zmin + z[k]*L;
		for (j=ys; j<ys+ym; j++) {
			Y = Ymin + y[j]*L;
			for (i=xs; i<xs+xm; i++) {
				X = Xmin + x[i]*L;
				if (bcType == 0 || bcType == 1 || bcType == 2 || bcType == 10 || bcType == 11) {
					// Defines upper and lower indices for the differenciation
					id.x[0] = i; id.x[1] = i-1; id.x[2] = i+1; 
					id.y[0] = j; id.y[1] = j-1; id.y[2] = j+1;
					id.z[0] = k; id.z[1] = k-1; id.z[2] = k+1;
					
					// Defines helpful coefficients for the differenciation
					if (i==0)         {dh.x[0] = -(x[id.x[2]] - x[id.x[0]]); dh.x[1] =   x[id.x[2]] - x[id.x[0]] ;}
					else if (i==mx-1) {dh.x[0] =   x[id.x[1]] - x[id.x[0]] ; dh.x[1] = -(x[id.x[1]] - x[id.x[0]]);}
					else              {dh.x[0] =   x[id.x[1]] - x[id.x[0]] ; dh.x[1] =   x[id.x[2]] - x[id.x[0]] ;}
					if (j==0)         {dh.y[0] = -(y[id.y[2]] - y[id.y[0]]); dh.y[1] =   y[id.y[2]] - y[id.y[0]] ;}
					else if (j==my-1) {dh.y[0] =   y[id.y[1]] - y[id.y[0]] ; dh.y[1] = -(y[id.y[1]] - y[id.y[0]]);}
					else              {dh.y[0] =   y[id.y[1]] - y[id.y[0]] ; dh.y[1] =   y[id.y[2]] - y[id.y[0]] ;}
					if (k==0)         {dh.z[0] = -(z[id.z[2]] - z[id.z[0]]); dh.z[1] =   z[id.z[2]] - z[id.z[0]] ;}
					else if (k==mz-1) {dh.z[0] =   z[id.z[1]] - z[id.z[0]] ; dh.z[1] = -(z[id.z[1]] - z[id.z[0]]);}
					else              {dh.z[0] =   z[id.z[1]] - z[id.z[0]] ; dh.z[1] =   z[id.z[2]] - z[id.z[0]] ;}
				}
				if (bcType == 3) {
					// Defines upper and lower indices for the differenciation
					if (i==0)         {id.x[0] = i; id.x[1] = i+1; id.x[2] = i+2;}
					else if (i==mx-1) {id.x[0] = i; id.x[1] = i-2; id.x[2] = i-1;}
					else              {id.x[0] = i; id.x[1] = i-1; id.x[2] = i+1;}
					if (j==0)         {id.y[0] = j; id.y[1] = j+1; id.y[2] = j+2;}
					else if (j==my-1) {id.y[0] = j; id.y[1] = j-2; id.y[2] = j-1;}
					else              {id.y[0] = j; id.y[1] = j-1; id.y[2] = j+1;}
					if (k==0)         {id.z[0] = k; id.z[1] = k+1; id.z[2] = k+2;}
					else if (k==mz-1) {id.z[0] = k; id.z[1] = k-2; id.z[2] = k-1;}
					else              {id.z[0] = k; id.z[1] = k-1; id.z[2] = k+1;}
					
					// Defines helpful coefficients for the differenciation
					dh.x[0] = x[id.x[1]] - x[id.x[0]]; dh.x[1] = x[id.x[2]] - x[id.x[0]];
					dh.y[0] = y[id.y[1]] - y[id.y[0]]; dh.y[1] = y[id.y[2]] - y[id.y[0]];
					dh.z[0] = z[id.z[1]] - z[id.z[0]]; dh.z[1] = z[id.z[2]] - z[id.z[0]];
				}
				
				// Calculate v on all interior points applies to all BC cases for v
				D1.x[0] = - 1.0/dh.x[0] - 1.0/dh.x[1]; 
				D1.x[1] =   1.0/dh.x[0] - 1.0/(dh.x[0]-dh.x[1]); 
				D1.x[2] =   1.0/dh.x[1] + 1.0/(dh.x[0]-dh.x[1]);
				D1.y[0] = - 1.0/dh.y[0] - 1.0/dh.y[1]; 
				D1.y[1] =   1.0/dh.y[0] - 1.0/(dh.y[0]-dh.y[1]); 
				D1.y[2] =   1.0/dh.y[1] + 1.0/(dh.y[0]-dh.y[1]);
				D1.z[0] = - 1.0/dh.z[0] - 1.0/dh.z[1]; 
				D1.z[1] =   1.0/dh.z[0] - 1.0/(dh.z[0]-dh.z[1]); 
				D1.z[2] =   1.0/dh.z[1] + 1.0/(dh.z[0]-dh.z[1]);
				D2.x[0] =   2.0/(dh.x[0] * dh.x[1]); 
				D2.x[1] =   2.0/(dh.x[0] * (dh.x[0]-dh.x[1])); 
				D2.x[2] =   2.0/(dh.x[1] * (dh.x[1]-dh.x[0]));
				D2.y[0] =   2.0/(dh.y[0] * dh.y[1]); 
				D2.y[1] =   2.0/(dh.y[0] * (dh.y[0]-dh.y[1])); 
				D2.y[2] =   2.0/(dh.y[1] * (dh.y[1]-dh.y[0]));
				D2.z[0] =   2.0/(dh.z[0] * dh.z[1]); 
				D2.z[1] =   2.0/(dh.z[0] * (dh.z[0]-dh.z[1])); 
				D2.z[2] =   2.0/(dh.z[1] * (dh.z[1]-dh.z[0]));
				
				// Calculate ne, ve, vi, J, E
				ni[O2p]     = u[k][j][i][d.ni[O2p]];     // nO2+
				ni[CO2p]    = u[k][j][i][d.ni[CO2p]];    // nCO2+
				ni[Op]      = u[k][j][i][d.ni[Op]];      // nO+
				ne = ni[O2p] + ni[CO2p] + ni[Op];
				vi[O2p][0]  = u[k][j][i][d.vi[O2p][0]];  // vxO2+
				vi[O2p][1]  = u[k][j][i][d.vi[O2p][1]];  // vyO2+
				vi[O2p][2]  = u[k][j][i][d.vi[O2p][2]];  // vzO2+
				vi[CO2p][0] = u[k][j][i][d.vi[CO2p][0]]; // vxCO2+
				vi[CO2p][1] = u[k][j][i][d.vi[CO2p][1]]; // vyCO2+
				vi[CO2p][2] = u[k][j][i][d.vi[CO2p][2]]; // vzCO2+
				vi[Op][0]   = u[k][j][i][d.vi[Op][0]];   // vxO+
				vi[Op][1]   = u[k][j][i][d.vi[Op][1]];   // vyO+
				vi[Op][2]   = u[k][j][i][d.vi[Op][2]];   // vzO+
				pe          = u[k][j][i][d.pe];          // pe
				
				B[0]        = u[k][j][i][d.B[0]];        // Bx
				B[1]        = u[k][j][i][d.B[1]];        // By
				B[2]        = u[k][j][i][d.B[2]];        // Bz
				
				// Intermediate calculations
				dpe.dx = D1.x[0]*u[k][j][id.x[0]][d.pe] + D1.x[1]*u[k][j][id.x[1]][d.pe] + D1.x[2]*u[k][j][id.x[2]][d.pe]; //dpe.dx 
				dpe.dy = D1.y[0]*u[k][id.y[0]][i][d.pe] + D1.y[1]*u[k][id.y[1]][i][d.pe] + D1.y[2]*u[k][id.y[2]][i][d.pe]; //dpe.dy 
				dpe.dz = D1.z[0]*u[id.z[0]][j][i][d.pe] + D1.z[1]*u[id.z[1]][j][i][d.pe] + D1.z[2]*u[id.z[2]][j][i][d.pe]; //dpe.dz 

				for (m=0; m<3; m++) {
					dB[m].dx = D1.x[0]*u[k][j][id.x[0]][d.B[m]] + D1.x[1]*u[k][j][id.x[1]][d.B[m]] + D1.x[2]*u[k][j][id.x[2]][d.B[m]]; //dB[m].dx 
					dB[m].dy = D1.y[0]*u[k][id.y[0]][i][d.B[m]] + D1.y[1]*u[k][id.y[1]][i][d.B[m]] + D1.y[2]*u[k][id.y[2]][i][d.B[m]]; //dB[m].dy 
					dB[m].dz = D1.z[0]*u[id.z[0]][j][i][d.B[m]] + D1.z[1]*u[id.z[1]][j][i][d.B[m]] + D1.z[2]*u[id.z[2]][j][i][d.B[m]]; //dB[m].dz
                    
                    			// check the values of these
				}
				
				// Current
				J[0] = dB[2].dy - dB[1].dz; 
				J[1] = dB[0].dz - dB[2].dx;
				J[2] = dB[1].dx - dB[0].dy;
				
				// Electron velocity
				ve[0] = ( (ni[O2p]*vi[O2p][0] + ni[CO2p]*vi[CO2p][0] + ni[Op]*vi[Op][0]) - J[0] )/ ne;
				ve[1] = ( (ni[O2p]*vi[O2p][1] + ni[CO2p]*vi[CO2p][1] + ni[Op]*vi[Op][1]) - J[1] )/ ne;
				ve[2] = ( (ni[O2p]*vi[O2p][2] + ni[CO2p]*vi[CO2p][2] + ni[Op]*vi[Op][2]) - J[2] )/ ne;

				// projectiles velocities
				for (m=0;m<3;m++) {
					vb[0][m] = un[m];    // CO2
					vb[1][m] = un[m];    // O
					vb[2][m] = vi[0][m]; // O2+
					vb[3][m] = vi[1][m]; // CO2+
					vb[4][m] = vi[2][m]; // O2+
					vb[5][m] = ve[m];    // e
				}
				
				// Temperature
				Te = Interpolate(user->RefProf, 8 ,Z, lin_flat);
				nn[CO2]     = Interpolate(user->RefProf, 5 ,Z, lin_flat);
				nn[O]       = Interpolate(user->RefProf, 4 ,Z, lin_flat);
				nu[CO2]     = Ven(CO2, nn[CO2], Te);
				nu[O]       = Ven(O  , nn[O]  , Te);

				Te = Interpolate(user->RefProf, 8 ,Z, lin_flat)/T0;
				nn[CO2]     = Interpolate(user->RefProf, 5 ,Z, lin_flat)/n0;
				nn[O]       = Interpolate(user->RefProf, 4 ,Z, lin_flat)/n0;
                
				nu[0]       = v41(nn[CO2] *n0 , Te      *T0);	// e    CO2
				nu[1]       = v42(nn[O]   *n0 , Te      *T0);	// e    O
				nu[2]       = v43(ni[O2p] *n0 , Te      *T0);	// e    O2+
				nu[3]       = v44(ni[CO2p]*n0 , Te      *T0);	// e    CO2+
				nu[4]       = v45(ni[Op]  *n0 , Te      *T0);	// e    O+
				nu[5]       = v46(ne      *n0 , Te      *T0);	// e    e
				
				// Elastic collisions term
				for (m=0;m<3;m++) {
					el_coll[m] = 0;
					for (b=0;b<6;b++) {
						el_coll[m] += nu[b]*tau*(vb[b][m]-ve[m]);
					}
				}

				// chemistry for Gen Ohms Law
				chem_nuS[e][0] = v1();                              // CO2 + hv
				chem_nuS[e][1] = 0;
				chem_nuS[e][2] = v3();                              // O + hv
				chem_nuS[e][3] = 2*v4(nn[O]     *n0 , Te    *T0);   // O + e
				for (m=0;m<3;m++) {
					E_chem_S[m]   = tau   *(chem_nuS[e][0]   *(un[m] - vi[e][m])
					+ chem_nuS[e][1]   *(un[m] - vi[e][m])
					+ chem_nuS[e][2]   *(un[m] - vi[e][m])
					+ chem_nuS[e][3]   *(un[m] - vi[e][m]));
				}
				
				// E-field
				if (vDamping) { //((1.0-tanh(Z-Lz)/(Lz/12.0))/2.0);
					for (m=0; m<3; m++) { 
						un[m]*=
						.5*(1+erf((X-XL)/lambda)) * .5*(1-erf((X-XU)/lambda)) *
						.5*(1+erf((Y-YL)/lambda)) * .5*(1-erf((Y-YU)/lambda)) *
						.5*(1+erf((Z-ZL)/lambda)) * .5*(1-erf((Z-ZU)/lambda)) ;
					}
				}

				// Gen Ohms Law
				E[0] = -CrossP(ve,B,0) - dpe.dx/ne;
				E[1] = -CrossP(ve,B,1) - dpe.dy/ne;
				E[2] = -CrossP(ve,B,2) - dpe.dz/ne;

				if(user->chemswitch==1){
					E[0] += E_chem_S[0];
					E[1] += E_chem_S[1];
					E[2] += E_chem_S[2];
				}
				if(user->collswitch==1){
					E[0] += el_coll[0];
					E[1] += el_coll[1];
					E[2] += el_coll[2];
				}
				
				// Store data
				for (m=0; m<3; m++) { 
					v[k][j][i][s.ve[m]] = ve[m];
					v[k][j][i][s.J[m]]  = J[m];
					v[k][j][i][s.E[m]]  = E[m];
				}
			}
		}
	}
	ierr = DMDAVecRestoreArrayDOF(db,V,&v);CHKERRQ(ierr);
	PetscFunctionReturn(0); 
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormFunction"
/* 
 FormFunction - Evaluates nonlinear function, F(u).
 
 Input Parameters:
 .  ts - the TS context
 .  U - input vector
 .  ctx - optional user-defined context, as set by SNESSetFunction()
 
 Output Parameter:
 .  F - function vector
 */
PetscErrorCode FormFunction(TS ts,PetscReal ftime,Vec U,Vec F,void *ctx)
{
	AppCtx         *user = (AppCtx*)ctx;
	PetscReal      n0 = user->n0, v0 = user->v0, p0 = user->p0, T0 = user->T0;
	DM             da = (DM)user->da;
	DM             db = (DM)user->db;
	PetscInt       mx = user->mx, my = user->my, mz = user->mz;
	PetscReal      L=user->L;
	PetscReal      Lz;
	PetscReal      *x=user->x,*y=user->y,*z=user->z;
	PetscInt       bcType = user->bcType;
	PetscBool      limiters = user->limiters;
	PetscBool      blockers = user->blockers;
	PetscBool      vDamping = user->vDamping;
	PetscReal      lambda=user->lambda;
	PetscReal      XL=user->inXmin,XU=user->inXmax;
	PetscReal      YL=user->inYmin,YU=user->inYmax;
	PetscReal      ZL=user->ZL    ,ZU=user->ZU    ;
	PetscReal      Xmin=user->outXmin,Ymin=user->outYmin, Zmin=user->outZmin;
	PetscReal      X,Y,Z;
	dvi            d = user->d;
	svi            s = user->s;
	PetscReal      rM = user->rM;
	PetscReal      tau = user->tau; 
	PetscReal      me = user->me;                                                  // me 
	PetscReal      mi[3] = {user->mi[0], user->mi[1], user->mi[2]};                // mi[O2+,CO2+,O+]
	PetscReal      un[3] = {user->un[0]/v0, user->un[1]/v0, user->un[2]/v0};       // neutral wind
	PetscReal      gama[4] = {user->gama[O2p],user->gama[CO2p],user->gama[Op],user->gama[e]};
	PetscReal      c=1.0/(PetscSqrtScalar(mu0*eps0)*v0);                           // Normalized speed of light (all points, all directions, all species
	PetscReal      vA[3] = {                                                       // Alfven speed
		0.0,                                                         // * 
		0.0,                                                         // * for a GIVEN species l
		0.0                                                          // * in ANY direction m
	};                                                                   // * at the CURRENT point (k,j,i)
	
	PetscReal      vA_max[3] = {                                                   // Maximum Alfven speed
		0.0,                                                         // * 
		0.0,                                                         // * for a GIVEN species l
		0.0                                                          // * in ANY direction m
	};                                                                   // * at the ENTIRE SUBDOMAIN
	
	PetscReal      vs[3] = {                                                       // Sonic speed
		0.0,                                                         // * 
		0.0,                                                         // * for a GIVEN species l
		0.0                                                          // * in ANY direction m
	};                                                                   // * at the CURRENT point (k,j,i)
	
	PetscReal      vs_max[3] = {                                                   // Maximum Sonic speed
		0.0,                                                         // * 
		0.0,                                                         // * for a GIVEN species l
		0.0                                                          // * in ANY direction m
	};                                                                   // * at the ENTIRE SUBDOMAIN
	
	PetscReal      vf[3] = {                                                       // Fast magnetosonic wave speed
		0.0,                                                         // * 
		0.0,                                                         // * for a GIVEN species l
		0.0                                                          // * in ANY direction m
	};                                                                   // * at the CURRENT point (k,j,i)
	
	PetscReal      vf_max[3] = {                                                   // Maximum Fast magnetosonic wave speed
		0.0,                                                         // * 
		0.0,                                                         // * for a GIVEN species l
		0.0                                                          // * in ANY direction m
	};                                                                   // * at the ENTIRE SUBDOMAIN
	
	PetscReal      ve_max[3] = {                                                   // Maximum speed of the hydrodynamic electron flow
		0.0,                                                         // * 
		0.0,                                                         // * 
		0.0                                                          // * in ANY direction m
	};                                                                   // * at the ENTIRE SUBDOMAIN
	
	PetscReal      vcr[3][3] = {                                                   // Critical speed
		{0.0, 0.0, 0.0},                                             // * sum of the hydrodynamic and fast magnetoacoustic wave speeds
		{0.0, 0.0, 0.0},                                             // * for a GIVEN species l
		{0.0, 0.0, 0.0}                                              // * in the GIVEN direction m
	};                                                                   // * at the CURRENT point (k,j,i)
	
	PetscReal      el_coll[3][3] = {                                               // Elastic collision term
		{0.0, 0.0, 0.0},                                             // * 
		{0.0, 0.0, 0.0},                                             // * for a GIVEN species l
		{0.0, 0.0, 0.0}                                              // * in the GIVEN direction m
	};                                                                   // * at the CURRENT point (k,j,i)
	
	PetscReal      vcr_max[3][3] = {                                               // Critical speed 
		{0.0, 0.0, 0.0},                                             // * sum of the hydro-dynamic and fast magnetoacoustic wave speeds
		{0.0, 0.0, 0.0},                                             // * for a GIVEN species l
		{0.0, 0.0, 0.0}                                              // * in the GIVEN direction m
	};                                                                   // * in the SUBDOMAIN (at the end of the space scanning of the subdomain)
	
	PetscReal      vi_cr[3][3] = {                                                 // max physical ion speed
		{0.0, 0.0, 0.0},                                             // *
		{0.0, 0.0, 0.0},                                             // * for a GIVEN species l
		{0.0, 0.0, 0.0}                                              // * in ANY direction m
	};
    
	PetscReal      vth = 0.0;                                                      // Thermal velocity
	
	PetscReal      T_cr[4] = {0.0,0.0,0.0,0.0};
	PetscErrorCode ierr;
	PetscReal      t,dt,dt_CFL,h[3];
	PetscReal      nu[4][6];                                                       // ion-neutral collision frequency
	PetscReal      cont_chem_S[4];
	PetscReal      cont_chem_L[4];
	PetscReal      mom_chem_S[3][3];
	PetscReal      ene_chem_S[4];
	PetscReal      ene_chem_L[4];
	PetscReal      chem_nuS[7][5];
	PetscReal      chem_nuL[7][5];
	PetscReal      gaman[2];                                                       // Specific heat of neutrals
	PetscReal      mn[2];                                                          // Mass of neutrals, kg
	PetscReal      E[3],B[3];                                                      // E-field, B-field
	PetscReal      ne,ni[3],nn[2];                                                 // [O2+], [CO2+], [O+], ne in m^-3
	PetscReal      Ti[3],Te,Tn[2];                                                 // Temperature of e, O2+, CO2+, O+, CO2, 0 in K
	PetscReal      pe,pi[3],pn[2];                                                       // pressure of O2+, CO2+, O+, e, O, and CO2 in J/m^3
	Vec            V;
	PetscReal      ****u,****v,****f;
	Vec            localU,localV;
	PetscInt       rank,numranks,step;
	PetscInt       b,i,j,k,l,m,n,xs,ys,zs,xm,ym,zm;
	PetscReal      ve[3],vi[3][3],vb[6][3];                                        // electron and ion velocities
	PetscReal      vdiffO2p[3],vdiffCO2p[3],vdiffOp[3],vdiffe[3];
	diff           dE[3],dpe,dve[3],dni[3],dpi[3],dvi[3][3];                       // Other derivatives 
	stencil        id;
	coeff3         D1;
	coeff2         dh;
	PetscBool      flag = PETSC_FALSE;
	FILE           *fp = NULL, *fp2 = NULL;
	char           fName[PETSC_MAX_PATH_LEN];
	char           fName2[PETSC_MAX_PATH_LEN];
	PetscInt       root = 0;
	PetscInt       tag;
	PetscInt       msgsource;
	PetscInt       kk;
    
	PetscFunctionBegin;
	id.x[0] = 0; id.y[0] = 0; id.z[0] = 0;
	id.x[1] = 0; id.y[1] = 0; id.z[1] = 0;
	id.x[2] = 0; id.y[2] = 0; id.z[2] = 0;
	
	dh.x[0] = 0; dh.y[0] = 0; dh.z[0] = 0;
	dh.x[1] = 0; dh.y[1] = 0; dh.z[1] = 0;
	
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &numranks);
	Lz = user->outZmax - user->outZmin;
	
	// Get current iteration t, dt, step
	ierr = TSGetTimeStep(ts,&dt);CHKERRQ(ierr);
	ierr = TSGetTime(ts,&t);CHKERRQ(ierr);
	ierr = TSGetStepNumber(ts,&step);CHKERRQ(ierr);
	
	ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
    
	ierr = DMGetGlobalVector(db,&V);CHKERRQ(ierr);
	
	// Create and fill ****u inner cells proc ghosts
	ierr = DMGetLocalVector(da,&localU);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(da,U,INSERT_VALUES,localU);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(da,U,INSERT_VALUES,localU);CHKERRQ(ierr);
	
	// Fill ****u domain ghosts
	ierr = DMDAVecGetArrayDOF(da,localU,&u);CHKERRQ(ierr);
	ierr = FormBCu(u,user);CHKERRQ(ierr);
	
	// Create and fill ****v inner cells
	ierr = FormIntermediateFunction(u,V,user);CHKERRQ(ierr);
	
	// Create and fill ****v inner cells proc ghosts
	ierr = DMGetLocalVector(db,&localV);CHKERRQ(ierr);
	
	ierr = DMGlobalToLocalBegin(db,V,INSERT_VALUES,localV);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(db,V,INSERT_VALUES,localV);CHKERRQ(ierr);
	
	// Fill ****v domain ghosts
	ierr = DMDAVecGetArrayDOF(db,localV,&v);CHKERRQ(ierr);
	ierr = FormBCv(v,user);CHKERRQ(ierr);
	
	// Check values sanity
	ierr = CheckValues(u,v,user);CHKERRQ(ierr);
	
	// Form F(u,v) = Udot
	ierr = DMDAVecGetArrayDOF(da,F,&f);CHKERRQ(ierr);

	if(step==0 || ceilf((step+1)/user->viz_dstep) == (float)(step+1)/user->viz_dstep){
		if(rank==0){
			// Commented by Kellen to remove .dat
			//sprintf(fName,"%s/Rates%d.dat",user->vName,step+1);
			//fp = fopen(fName,"w");
		}
	}

	// Compute function over the locally owned part of the grid

	dt_CFL = 1e300;
	flag   = PETSC_FALSE;
	for (k=zs; k<zs+zm; k++) {
		Z = Zmin + z[k]*L;
		for (j=ys; j<ys+ym; j++) {
			Y = Ymin + y[j]*L;
			for (i=xs; i<xs+xm; i++) {
				X = Xmin + x[i]*L;
				if (bcType == 0 || bcType == 1 || bcType == 2 || bcType == 10 || bcType == 11) {
					// Defines upper and lower indices for the differenciation
					id.x[0] = i; id.x[1] = i-1; id.x[2] = i+1; 
					id.y[0] = j; id.y[1] = j-1; id.y[2] = j+1;
					id.z[0] = k; id.z[1] = k-1; id.z[2] = k+1;
					
					// Defines helpful coefficients for the differenciation
					if (i==0)         {dh.x[0] = -(x[id.x[2]] - x[id.x[0]]); dh.x[1] =   x[id.x[2]] - x[id.x[0]] ;}
					else if (i==mx-1) {dh.x[0] =   x[id.x[1]] - x[id.x[0]] ; dh.x[1] = -(x[id.x[1]] - x[id.x[0]]);}
					else              {dh.x[0] =   x[id.x[1]] - x[id.x[0]] ; dh.x[1] =   x[id.x[2]] - x[id.x[0]] ;}
					if (j==0)         {dh.y[0] = -(y[id.y[2]] - y[id.y[0]]); dh.y[1] =   y[id.y[2]] - y[id.y[0]] ;}
					else if (j==my-1) {dh.y[0] =   y[id.y[1]] - y[id.y[0]] ; dh.y[1] = -(y[id.y[1]] - y[id.y[0]]);}
					else              {dh.y[0] =   y[id.y[1]] - y[id.y[0]] ; dh.y[1] =   y[id.y[2]] - y[id.y[0]] ;}
					if (k==0)         {dh.z[0] = -(z[id.z[2]] - z[id.z[0]]); dh.z[1] =   z[id.z[2]] - z[id.z[0]] ;}
					else if (k==mz-1) {dh.z[0] =   z[id.z[1]] - z[id.z[0]] ; dh.z[1] = -(z[id.z[1]] - z[id.z[0]]);}
					else              {dh.z[0] =   z[id.z[1]] - z[id.z[0]] ; dh.z[1] =   z[id.z[2]] - z[id.z[0]] ;}
				}
				if (bcType == 3) {
					// Defines upper and lower indices for the differenciation
					if (i==0)         {id.x[0] = i; id.x[1] = i+1; id.x[2] = i+2;}
					else if (i==mx-1) {id.x[0] = i; id.x[1] = i-2; id.x[2] = i-1;}
					else              {id.x[0] = i; id.x[1] = i-1; id.x[2] = i+1;}
					if (j==0)         {id.y[0] = j; id.y[1] = j+1; id.y[2] = j+2;}
					else if (j==my-1) {id.y[0] = j; id.y[1] = j-2; id.y[2] = j-1;}
					else              {id.y[0] = j; id.y[1] = j-1; id.y[2] = j+1;}
					if (k==0)         {id.z[0] = k; id.z[1] = k+1; id.z[2] = k+2;}
					else if (k==mz-1) {id.z[0] = k; id.z[1] = k-2; id.z[2] = k-1;}
					else              {id.z[0] = k; id.z[1] = k-1; id.z[2] = k+1;}
					
					// Defines helpful coefficients for the differenciation
					dh.x[0] = x[id.x[1]] - x[id.x[0]]; dh.x[1] = x[id.x[2]] - x[id.x[0]];
					dh.y[0] = y[id.y[1]] - y[id.y[0]]; dh.y[1] = y[id.y[2]] - y[id.y[0]];
					dh.z[0] = z[id.z[1]] - z[id.z[0]]; dh.z[1] = z[id.z[2]] - z[id.z[0]];
				}
                
				
				// Coefficients for differenciation
				D1.x[0]  = - 1.0/dh.x[0] - 1.0/dh.x[1]; 
				D1.x[1]  =   1.0/dh.x[0] - 1.0/(dh.x[0]-dh.x[1]); 
				D1.x[2]  =   1.0/dh.x[1] + 1.0/(dh.x[0]-dh.x[1]);
				D1.y[0]  = - 1.0/dh.y[0] - 1.0/dh.y[1]; 
				D1.y[1]  =   1.0/dh.y[0] - 1.0/(dh.y[0]-dh.y[1]); 
				D1.y[2]  =   1.0/dh.y[1] + 1.0/(dh.y[0]-dh.y[1]);
				D1.z[0]  = - 1.0/dh.z[0] - 1.0/dh.z[1]; 
				D1.z[1]  =   1.0/dh.z[0] - 1.0/(dh.z[0]-dh.z[1]); 
				D1.z[2]  =   1.0/dh.z[1] + 1.0/(dh.z[0]-dh.z[1]);
				
				// For readability, we distribute the components of u onto the variable n, v, p, B
				ni[O2p]     = u[k][j][i][d.ni[O2p]];     // nO2+
				ni[CO2p]    = u[k][j][i][d.ni[CO2p]];    // nCO2+
				ni[Op]      = u[k][j][i][d.ni[Op]];      // nO+
				ne = ni[O2p] + ni[CO2p] + ni[Op];
				vi[O2p][0]  = u[k][j][i][d.vi[O2p][0]];  // vxO2+
				vi[O2p][1]  = u[k][j][i][d.vi[O2p][1]];  // vyO2+
				vi[O2p][2]  = u[k][j][i][d.vi[O2p][2]];  // vzO2+
				vi[CO2p][0] = u[k][j][i][d.vi[CO2p][0]]; // vxCO2+
				vi[CO2p][1] = u[k][j][i][d.vi[CO2p][1]]; // vyCO2+
				vi[CO2p][2] = u[k][j][i][d.vi[CO2p][2]]; // vzCO2+
				vi[Op][0]   = u[k][j][i][d.vi[Op][0]];   // vxO+
				vi[Op][1]   = u[k][j][i][d.vi[Op][1]];   // vyO+
				vi[Op][2]   = u[k][j][i][d.vi[Op][2]];   // vzO+
				pi[O2p]     = u[k][j][i][d.pi[O2p]];     // pO2+
				pi[CO2p]    = u[k][j][i][d.pi[CO2p]];    // pCO2+
				pi[Op]      = u[k][j][i][d.pi[Op]];      // pO+
				pe          = u[k][j][i][d.pe];          // pe
				
				B[0]        = u[k][j][i][d.B[0]];        // Bx
				B[1]        = u[k][j][i][d.B[1]];        // By
				B[2]        = u[k][j][i][d.B[2]];        // Bz
				
				ve[0]       = v[k][j][i][s.ve[0]];       // vxe
				ve[1]       = v[k][j][i][s.ve[1]];       // vye
				ve[2]       = v[k][j][i][s.ve[2]];       // vze
				E[0]        = v[k][j][i][s.E[0]];        // Ex
				E[1]        = v[k][j][i][s.E[1]];        // Ey
				E[2]        = v[k][j][i][s.E[2]];        // Ez
				
                
				nn[CO2] = Interpolate(user->RefProf, 5,Z,lin_exp )/n0;
				nn[O]   = Interpolate(user->RefProf, 4,Z,lin_exp )/n0;
				for (l=0; l<3; l++)
					Ti[l] = Interpolate(user->RefProf, 9 ,Z, lin_flat)/T0;
				Tn[CO2] = Interpolate(user->RefProf,10,Z,lin_flat)/T0;
				Tn[O]   = Interpolate(user->RefProf,10,Z,lin_flat)/T0;
				Te = Interpolate(user->RefProf, 8 ,Z, lin_flat)/T0;
                
				// projectiles velocities
				for (m=0;m<3;m++) {
					vb[0][m] = un[m];    // CO2
					vb[1][m] = un[m];    // O
					vb[2][m] = vi[0][m]; // O2+
					vb[3][m] = vi[1][m]; // CO2+
					vb[4][m] = vi[2][m]; // O2+
					vb[5][m] = ve[m];    // e
				}
				
				// collisions
				nu[0][0] = v11(nn[CO2] *n0);                              // O2+  CO2
				nu[0][1] = v12(nn[O]   *n0);                              // O2+  O
				nu[0][2] = v13(ni[O2p] *n0 , Ti[O2p] *T0);                // O2+  O2+
				nu[0][3] = v14(ni[CO2p]*n0 , Ti[CO2p]*T0);                // O2+  CO2+
				nu[0][4] = v15(ni[Op]  *n0 , Ti[Op]  *T0);                // O2+  O+
				nu[0][5] = v16();                                         // O2+  e
				
				nu[1][0] = v21(nn[CO2] *n0 , Tn[CO2] *T0 , Ti[CO2p] *T0); // CO2+ CO2
				nu[1][1] = v22(nn[O]   *n0);                              // CO2+ O
				nu[1][2] = v23(ni[O2p] *n0 , Ti[O2p] *T0);                // CO2+ O2+
				nu[1][3] = v24(ni[CO2p]*n0 , Ti[CO2p]*T0);                // CO2+ CO2+
				nu[1][4] = v25(ni[Op]  *n0 , Ti[Op]  *T0);                // CO2+ O+
				nu[1][5] = v26();                                         // CO2+ e
				
				nu[2][0] = v31(nn[CO2] *n0);                              // O+   CO2
				nu[2][1] = v32(nn[O]   *n0 , Tn[O]   *T0 , Ti[Op]   *T0); // O+   O
				nu[2][2] = v33(ni[O2p] *n0 , Ti[O2p] *T0);                // O+   O2+
				nu[2][3] = v34(ni[CO2p]*n0 , Ti[CO2p]*T0);                // O+   CO2+
				nu[2][4] = v35(ni[Op]  *n0 , Ti[Op]  *T0);                // O+   O+
				nu[2][5] = v36();                                         // O+   e
				
				nu[3][0] = v41(nn[CO2] *n0 , Te      *T0);                // e    CO2
				nu[3][1] = v42(nn[O]   *n0 , Te      *T0);                // e    O
				nu[3][2] = v43(ni[O2p] *n0 , Te      *T0);                // e    O2+
				nu[3][3] = v44(ni[CO2p]*n0 , Te      *T0);                // e    CO2+
				nu[3][4] = v45(ni[Op]  *n0 , Te      *T0);                // e    O+
				nu[3][5] = v46(ne      *n0 , Te      *T0);                // e    e
                
				for (m=0;m<3;m++) {
					for (l=0;l<3;l++) {
						el_coll[l][m] = 0;
						for (b=0;b<6;b++) {
							el_coll[l][m] += tau*nu[l][b]*(vb[b][m]-vi[l][m]);  // Elastic Collisions for momentum,
						}
					}
				}
                
 				// CHEMISTRY RATE CALCULATIONS
				chem_nuS[O2p][0] = v8(ni[CO2p]  *n0);               // CO2+ + O
				chem_nuS[O2p][1] = v10(ni[Op]   *n0);               // O+ + CO2
				chem_nuL[O2p][0] = v5(ne        *n0 , Te    *T0);   // O2+ + e
                
				chem_nuS[CO2p][0] = v1();                           // CO2 + hv
				chem_nuS[CO2p][1] = 0;
				chem_nuL[CO2p][0] = v6(ne       *n0 , Te    *T0);   // CO2+ + e
				chem_nuL[CO2p][1] = v8(nn[O]    *n0);               // CO2+ + O
				chem_nuL[CO2p][2] = v9(nn[O]    *n0);               // CO2+ + O
                
				chem_nuS[Op][0] = v3();                             // O + hv
				chem_nuS[Op][1] = v4(ne         *n0 , Te    *T0);   // O + e
				chem_nuS[Op][2] = v9(ni[CO2p]   *n0);               // O + CO2+
				chem_nuL[Op][0] = v7(ne         *n0 , Te    *T0);   // O+ + e
				chem_nuL[Op][1] = v10(nn[CO2]   *n0);               // O+ + CO2
                
				chem_nuS[e][0] = v1();                              // CO2 + hv
				chem_nuS[e][1] = 0;
				chem_nuS[e][2] = v3();                              // O + hv
				chem_nuS[e][3] = 2*v4(nn[O]     *n0 , Te    *T0);   // O + e
				chem_nuL[e][0] = 0;
				chem_nuL[e][1] = v4(nn[O]       *n0 , Te    *T0);   // O + e
				chem_nuL[e][2] = v5(ni[O2p]     *n0 , Te    *T0);   // O2+ + e
				chem_nuL[e][3] = v6(ni[CO2p]    *n0 , Te    *T0);   // CO2+ + e
				chem_nuL[e][4] = v7(ni[Op]      *n0 , Te    *T0);   // O+ + e
            
				// CONSERVATION OF MASS SOURCES AND LOSSES
				cont_chem_S[O2p]  = tau     *(chem_nuS[O2p][0] *nn[O]    + chem_nuS[O2p][1] *nn[CO2]);
				cont_chem_S[CO2p] = tau     *(chem_nuS[CO2p][0]*nn[CO2]  + chem_nuS[CO2p][1]*nn[CO2]);
				cont_chem_S[Op]   = tau     *(chem_nuS[Op][0]  *nn[O]    + chem_nuS[Op][1]  *nn[O]     + chem_nuS[Op][2]  *nn[O]);
				cont_chem_L[O2p]  = tau     *(chem_nuL[O2p][0] *ni[Op]);
				cont_chem_L[CO2p] = tau     *(chem_nuL[CO2p][0]*ni[CO2p] + chem_nuL[CO2p][1]*ni[CO2p]  + chem_nuL[CO2p][2]*ni[CO2p]);
				cont_chem_L[Op]   = tau     *(chem_nuL[Op][0]  *ni[Op]   + chem_nuL[Op][1]  *ni[Op]);
                
				// CONSERVATION OF MOMENTUM SOURCES
				for (m=0;m<3;m++) {
					mom_chem_S[O2p][m]  = tau   *(chem_nuS[O2p][0]  *(nn[O]/ni[O2p])    *(un[m] - vi[O2p][m])
						+ chem_nuS[O2p][1]  *(nn[CO2]/ni[O2p])  *(un[m] - vi[O2p][m]));
					mom_chem_S[CO2p][m] = tau   *(chem_nuS[CO2p][0] *(nn[CO2]/ni[CO2p]) *(un[m] - vi[CO2p][m])
						+ chem_nuS[CO2p][1] *(nn[CO2]/ni[CO2p]) *(un[m] - vi[CO2p][m]));
					mom_chem_S[Op][m]   = tau   *(chem_nuS[Op][0]   *(nn[O]/ni[Op])     *(un[m] - vi[Op][m])
						+ chem_nuS[Op][1]   *(nn[O]/ni[Op])     *(un[m] - vi[Op][m])
						+ chem_nuS[Op][2]   *(nn[O]/ni[Op])     *(un[m] - vi[Op][m]));
				}
                
				// CONSERVATION OF ENERGY SOURCES AND LOSSES
				mn[CO2]      = 1.9944235E-26 + 2*2.6566962E-26;
				mn[O]        = 2.6566962E-26;
				pn[CO2]      = nn[CO2]*n0*kB*Tn[CO2]*T0/p0;
				pn[O]        = nn[O]*n0*kB*Tn[O]*T0/p0;
				gaman[CO2]   = 1.289;
				gaman[O]     = 1.667;
				for (m=0;m<3;m++) {
					vdiffO2p[m]  = un[m] - vi[O2p][m];
					vdiffCO2p[m] = un[m] - vi[CO2p][m];
					vdiffOp[m]   = un[m] - vi[Op][m];
					vdiffe[m]    = un[m] - ve[m];
				}
				ene_chem_S[O2p]     = tau   *(chem_nuS[O2p][0]      *(mi[O2p]/mn[O])    *(gama[O2p]-1)/(gaman[O]-1)     *pn[O]
					+ chem_nuS[O2p][1]      *(mi[O2p]/mn[CO2])  *(gama[O2p]-1)/(gaman[CO2]-1)   *pn[CO2]
					+ (gama[O2p]-1)   *chem_nuS[O2p][0]    *(mi[O2p]*nn[O])/me    *Norm2(vdiffCO2p)/2
					+ (gama[O2p]-1)   *chem_nuS[O2p][1]    *(mi[O2p]*nn[CO2])/me  *Norm2(vdiffCO2p)/2);
				ene_chem_S[CO2p]    = tau   *(chem_nuS[CO2p][0]     *(mi[CO2p]/mn[CO2]) *(gama[CO2p]-1)/(gaman[O]-1)    *pn[CO2]
					+ chem_nuS[CO2p][1]     *(mi[CO2p]/mn[CO2]) *(gama[CO2p]-1)/(gaman[CO2]-1)  *pn[CO2]
					+ (gama[CO2p]-1)  *chem_nuS[CO2p][0]   *(mi[CO2p]*nn[CO2])/me *Norm2(vdiffCO2p)/2
					+ (gama[CO2p]-1)  *chem_nuS[CO2p][1]   *(mi[CO2p]*nn[CO2])/me *Norm2(vdiffCO2p)/2);
				ene_chem_S[Op]      = tau   *(chem_nuS[Op][0]       *(mi[Op]/mn[O])     *(gama[Op]-1)/(gaman[O]-1)      *pn[O]
					+ chem_nuS[Op][1]       *(mi[Op]/mn[O])     *(gama[Op]-1)/(gaman[O]-1)    *pn[O]
					+ chem_nuS[Op][2]       *(mi[Op]/mn[O])     *(gama[Op]-1)/(gaman[O]-1)    *pn[O]
					+ (gama[Op]-1)    *chem_nuS[Op][0]     *(mi[Op]*nn[O])/me     *Norm2(vdiffOp)/2
					+ (gama[Op]-1)    *chem_nuS[Op][1]     *(mi[Op]*nn[O])/me     *Norm2(vdiffOp)/2
					+ (gama[Op]-1)    *chem_nuS[Op][2]     *(mi[Op]*nn[O])/me     *Norm2(vdiffOp)/2);
				ene_chem_S[e]       = tau   *(chem_nuS[e][0]       *(mi[Op]/mn[CO2])   *(gama[e]-1)/(gaman[CO2]-1)    *pn[O]
					+ chem_nuS[e][1]       *(mi[Op]/mn[CO2])   *(gama[e]-1)/(gaman[CO2]-1)    *pn[O]
					+ chem_nuS[e][2]       *(mi[Op]/mn[O])     *(gama[e]-1)/(gaman[O]-1)      *pn[O]
					+ chem_nuS[e][3]       *(mi[Op]/mn[O])     *(gama[e]-1)/(gaman[O]-1)      *pn[O]
					+ (gama[e]-1)    *chem_nuS[e][0]     *(me*nn[CO2])/me         *Norm2(vdiffe)/2
					+ (gama[e]-1)    *chem_nuS[e][1]     *(me*nn[CO2])/me         *Norm2(vdiffe)/2
					+ (gama[e]-1)    *chem_nuS[e][2]     *(me*nn[O])/me           *Norm2(vdiffe)/2
					+ (gama[e]-1)    *chem_nuS[e][3]     *(me*nn[O])/me           *Norm2(vdiffe)/2);
				ene_chem_L[O2p]     = tau   *(chem_nuL[O2p][0] *pi[Op]);
				ene_chem_L[CO2p]    = tau   *(chem_nuL[CO2p][0]*pi[CO2p] + chem_nuL[CO2p][1]*pi[CO2p]  + chem_nuL[CO2p][2]*pi[CO2p]);
				ene_chem_L[Op]      = tau   *(chem_nuL[Op][0]  *pi[Op]   + chem_nuL[Op][1]  *pi[Op]);
				ene_chem_L[e]       = tau   *(chem_nuL[e][0]   *pe       + chem_nuL[e][1]   *pe        + chem_nuL[e][2]   *pe
					+ chem_nuL[e][3] *pe       + chem_nuL[e][4]   *pe);
				
				// Intermediate calculations
				for (m=0; m<3; m++) {
					dE[m].dx = D1.x[0]*v[k][j][id.x[0]][s.E[m]] + D1.x[1]*v[k][j][id.x[1]][s.E[m]] + D1.x[2]*v[k][j][id.x[2]][s.E[m]]; //dE[m].dx 
					dE[m].dy = D1.y[0]*v[k][id.y[0]][i][s.E[m]] + D1.y[1]*v[k][id.y[1]][i][s.E[m]] + D1.y[2]*v[k][id.y[2]][i][s.E[m]]; //dE[m].dy 
					dE[m].dz = D1.z[0]*v[id.z[0]][j][i][s.E[m]] + D1.z[1]*v[id.z[1]][j][i][s.E[m]] + D1.z[2]*v[id.z[2]][j][i][s.E[m]]; //dE[m].dz 
					
					// Do the same thing here for B as for E
                    
					dve[m].dx = D1.x[0]*v[k][j][id.x[0]][s.ve[m]] + D1.x[1]*v[k][j][id.x[1]][s.ve[m]] + D1.x[2]*v[k][j][id.x[2]][s.ve[m]]; //dve[m].dx 
					dve[m].dy = D1.y[0]*v[k][id.y[0]][i][s.ve[m]] + D1.y[1]*v[k][id.y[1]][i][s.ve[m]] + D1.y[2]*v[k][id.y[2]][i][s.ve[m]]; //dve[m].dy 
					dve[m].dz = D1.z[0]*v[id.z[0]][j][i][s.ve[m]] + D1.z[1]*v[id.z[1]][j][i][s.ve[m]] + D1.z[2]*v[id.z[2]][j][i][s.ve[m]]; //dve[m].dz 
				}
				dpe.dx = D1.x[0]*u[k][j][id.x[0]][d.pe] + D1.x[1]*u[k][j][id.x[1]][d.pe] + D1.x[2]*u[k][j][id.x[2]][d.pe]; //dpe.dx 
				dpe.dy = D1.y[0]*u[k][id.y[0]][i][d.pe] + D1.y[1]*u[k][id.y[1]][i][d.pe] + D1.y[2]*u[k][id.y[2]][i][d.pe]; //dpe.dy 
				dpe.dz = D1.z[0]*u[id.z[0]][j][i][d.pe] + D1.z[1]*u[id.z[1]][j][i][d.pe] + D1.z[2]*u[id.z[2]][j][i][d.pe]; //dpe.dz 
				
				for (l=0; l<3; l++) {
					dni[l].dx = D1.x[0]*u[k][j][id.x[0]][d.ni[l]] + D1.x[1]*u[k][j][id.x[1]][d.ni[l]] + D1.x[2]*u[k][j][id.x[2]][d.ni[l]]; //dni.dx 
					dni[l].dy = D1.y[0]*u[k][id.y[0]][i][d.ni[l]] + D1.y[1]*u[k][id.y[1]][i][d.ni[l]] + D1.y[2]*u[k][id.y[2]][i][d.ni[l]]; //dni.dy 
					dni[l].dz = D1.z[0]*u[id.z[0]][j][i][d.ni[l]] + D1.z[1]*u[id.z[1]][j][i][d.ni[l]] + D1.z[2]*u[id.z[2]][j][i][d.ni[l]]; //dni.dz 
					
					for (m=0; m<3; m++) {
						dvi[l][m].dx = D1.x[0]*u[k][j][id.x[0]][d.vi[l][m]] + D1.x[1]*u[k][j][id.x[1]][d.vi[l][m]] + D1.x[2]*u[k][j][id.x[2]][d.vi[l][m]]; //dvi[m].dx 
						dvi[l][m].dy = D1.y[0]*u[k][id.y[0]][i][d.vi[l][m]] + D1.y[1]*u[k][id.y[1]][i][d.vi[l][m]] + D1.y[2]*u[k][id.y[2]][i][d.vi[l][m]]; //dvi[m].dy 
						dvi[l][m].dz = D1.z[0]*u[id.z[0]][j][i][d.vi[l][m]] + D1.z[1]*u[id.z[1]][j][i][d.vi[l][m]] + D1.z[2]*u[id.z[2]][j][i][d.vi[l][m]]; //dvi[m].dz 
					}
					
					dpi[l].dx = D1.x[0]*u[k][j][id.x[0]][d.pi[l]] + D1.x[1]*u[k][j][id.x[1]][d.pi[l]] + D1.x[2]*u[k][j][id.x[2]][d.pi[l]]; //dpi.dx 
					dpi[l].dy = D1.y[0]*u[k][id.y[0]][i][d.pi[l]] + D1.y[1]*u[k][id.y[1]][i][d.pi[l]] + D1.y[2]*u[k][id.y[2]][i][d.pi[l]]; //dpi.dy 
					dpi[l].dz = D1.z[0]*u[id.z[0]][j][i][d.pi[l]] + D1.z[1]*u[id.z[1]][j][i][d.pi[l]] + D1.z[2]*u[id.z[2]][j][i][d.pi[l]]; //dpi.dz 
				}
				
				// CFL
				h[0] = MinAbs(dh.x[0],dh.x[1]);
				h[1] = MinAbs(dh.y[0],dh.y[1]);
				h[2] = MinAbs(dh.z[0],dh.z[1]);
				
				// Critical speed and time for ions
				// group v and phase v, wave is no physical? ask nykyri or ma
				for (l=0; l<3; l++) {
					vs[l]       = PetscSqrtScalar( (gama[e]*pe + gama[l]*pi[l]) / (ni[l]*mi[l]/me) );   // Sonic speed
					vs_max[l]   = MaxAbs( vs_max[l], vs[l] );                                           // Maximum Sonic speed
					
					vA[l]       = PetscSqrtScalar( Norm2(B) / (ni[l]*mi[l]/me) );                       // Alfven speed
					vA_max[l]   = MaxAbs( vA_max[l], vA[l] );                                           // Maximum Alfven speed
					
					vf[l]       = c*PetscSqrtScalar( (vs[l]*vs[l]+vA[l]*vA[l]) / (c*c+vA[l]*vA[l]) );   // Fast magnetosonic wave speed
					vf_max[l]   = MaxAbs( vf_max[l], vf[l] );                                           // Maximum Fast magnetosonic wave speed
					
					for (m=0; m<3; m++) {
						vcr_max[l][m]   = MaxAbs( vcr_max[l][m], PetscAbsScalar(vi[l][m]) );   // Max critical speed
					}

					/* T_cr[l]:
					 * critical time
					 * in all the subdomain up to the current point (k,j,i)
					 * for each species l
					 * indepently from all the direction in space (usually denoted m)
					 */

					T_cr[l] = MinAbs(MinAbs(h[0]/vcr_max[l][0] , h[1]/vcr_max[l][1]),h[2]/vcr_max[l][2]);
				}
				
				vth = sqrt(3.0/2.0*kB*Te/me)/v0;

				// Critical speed and time for electrons
				for (m=0; m<3; m++) {
					ve_max[m] = MaxAbs(MaxAbs(ve_max[m], ve[m]),vth);
					for (l=0; l<3; l++) {
						if (ve_max[m] > vcr_max[l][m]) { 
							//Case when the electron related velocities are the critical velocities
							if(blockers) {
								/* 
								 * If flag blockers is TRUE:
								 * End the run
								 * Normally ve>vi
								 */
								PetscPrintf(PETSC_COMM_WORLD,"ve_max[%.3d %.3d %.3d]=%12.6e\nvth[%.3d %.3d %.3d]=%12.6e\nvcr_max[%.3d %.3d %.3d][%d][%d]=%12.6e\n",k,j,i,ve_max[m]*v0,k,j,i,vth*v0,k,j,i,l,m,vcr_max[l][m]*v0);
								PetscPrintf(PETSC_COMM_WORLD,"The CFL criterion is not stringent enough!\n");
								exit(12);
							} else if (limiters) {
								/*
								 * If flag limiters is TRUE
								 * Adjust the critical time to account for the electron velocities
								 */
								T_cr[e] = MinAbs(MinAbs(h[0]/ve_max[0] , h[1]/ve_max[1]),h[2]/ve_max[2]);
							} else {
								/*
								 * Limiters and Blockers flags are FALSE
								 * Ignore electron velocities in the definition of the critical time for CFL criterion
								 * Set arbitrarily large value for T_cr[e]
								 */
								T_cr[e] = 1.0e+308;
							}
						}
					}
				}
                
				/* dt_CFL:
				 * critical time for CFL condition
				 * in all the subdomain up to the current point (k,j,i)
				 * for ALL species (usually denoted l)
				 * indepently from ANY the direction in space (usually denoted m)
				 */
				dt_CFL = MinAbs(MinAbs(MinAbs(T_cr[O2p],T_cr[CO2p]),T_cr[Op]),T_cr[e]);
				
				// Equations
				// Eq 0-2: Continuity equations for ions //
				for (l=0; l<3; l++){
					f[k][j][i][d.ni[l]] = -ni[l]*(dvi[l][0].dx + dvi[l][1].dy + dvi[l][2].dz) - (vi[l][0]*dni[l].dx + vi[l][1]*dni[l].dy + vi[l][2]*dni[l].dz);
					if(user->chemswitch==1){
						f[k][j][i][d.ni[l]] += cont_chem_S[l] - cont_chem_L[l];
					}
				}
        
				// Eq 3-11: Ion momenta // NEGLECT O and e-n collisions FOR NOW (07-29-11)
				if (vDamping) { //((1.0-tanh(Z-Lz)/(Lz/12.0))/2.0);
					for (m=0; m<3; m++) { 
						un[m]*=
						.5*(1+erf((X-XL)/lambda)) * .5*(1-erf((X-XU)/lambda)) *
						.5*(1+erf((Y-YL)/lambda)) * .5*(1-erf((Y-YU)/lambda)) *
						.5*(1+erf((Z-ZL)/lambda)) * .5*(1-erf((Z-ZU)/lambda)) ;
					}
				}
                
				// Eq 3-11: Ion momenta //
				for (l=0; l<3; l++) {
					f[k][j][i][d.vi[l][0]] = - (vi[l][0]*dvi[l][0].dx +vi[l][1]*dvi[l][0].dy +vi[l][2]*dvi[l][0].dz) + me/mi[l]*(E[0] +CrossP(vi[l],B,0) -dpi[l].dx/ni[l]);
					f[k][j][i][d.vi[l][1]] = - (vi[l][0]*dvi[l][1].dx +vi[l][1]*dvi[l][1].dy +vi[l][2]*dvi[l][1].dz) + me/mi[l]*(E[1] +CrossP(vi[l],B,1) -dpi[l].dy/ni[l]);
					f[k][j][i][d.vi[l][2]] = - (vi[l][0]*dvi[l][2].dx +vi[l][1]*dvi[l][2].dy +vi[l][2]*dvi[l][2].dz) + me/mi[l]*(E[2] +CrossP(vi[l],B,2) -dpi[l].dz/ni[l]);
					if(user->gravswitch==1)
						f[k][j][i][d.vi[l][2]] += - PetscPowScalar(1+Z/rM,-2.0);
					if(user->chemswitch==1){
						f[k][j][i][d.vi[l][0]] += mom_chem_S[l][0];
						f[k][j][i][d.vi[l][1]] += mom_chem_S[l][1];
						f[k][j][i][d.vi[l][2]] += mom_chem_S[l][2];
					}
					if(user->collswitch==1){
						f[k][j][i][d.vi[l][0]] += el_coll[l][0];
						f[k][j][i][d.vi[l][1]] += el_coll[l][1];
						f[k][j][i][d.vi[l][2]] += el_coll[l][2];
					}
				}

				// Eq 12-14: Equations of state for ions
				// Read more about heat transfer (Schunk and Nagy Book). Appendix I. page 518
				for (l=0; l<3; l++){
					f[k][j][i][d.pi[l]] = - gama[l] * pi[l] * (dvi[l][0].dx + dvi[l][1].dy + dvi[l][2].dz) - (vi[l][0]*dpi[l].dx + vi[l][1]*dpi[l].dy + vi[l][2]*dpi[l].dz);
					if(user->chemswitch==1){
						f[k][j][i][d.pi[l]] += ene_chem_S[l] - ene_chem_L[l];
					}
				}
				
				// Eq 15: Equation of state for electrons
				f[k][j][i][d.pe] = - gama[l] * pe * (dve[0].dx + dve[1].dy + dve[2].dz ) - ( ve[0]*dpe.dx + ve[1]*dpe.dy + ve[2]*dpe.dz);

				if(user->chemswitch==1){
					f[k][j][i][d.pe] += ene_chem_S[e] - ene_chem_L[e];
				}
				
				// Eq 16-18: Faraday's law
				// Ensure div B = 0 needs to be modded here
				f[k][j][i][d.B[0]] = - (dE[2].dy - dE[1].dz); 
				f[k][j][i][d.B[1]] = - (dE[0].dz - dE[2].dx);
				f[k][j][i][d.B[2]] = - (dE[1].dx - dE[0].dy);
			}
		}
	}

	m = MPI_Allreduce(MPI_IN_PLACE,&vs_max[0]  ,3 ,MPIU_REAL,MPIU_MAX,PETSC_COMM_WORLD);
	m = MPI_Allreduce(MPI_IN_PLACE,&vA_max[0]  ,3 ,MPIU_REAL,MPIU_MAX,PETSC_COMM_WORLD);
	m = MPI_Allreduce(MPI_IN_PLACE,&vf_max[0]  ,3 ,MPIU_REAL,MPIU_MAX,PETSC_COMM_WORLD);
	m = MPI_Allreduce(MPI_IN_PLACE,&ve_max[0]  ,3 ,MPIU_REAL,MPIU_MAX,PETSC_COMM_WORLD);
	m = MPI_Allreduce(MPI_IN_PLACE,&vcr_max[0] ,9 ,MPIU_REAL,MPIU_MAX,PETSC_COMM_WORLD);
	m = MPI_Allreduce(MPI_IN_PLACE,&dt_CFL ,1 ,MPIU_REAL,MPIU_MIN,PETSC_COMM_WORLD);
	
	// Adapt timestep to satisfy CFL
	user->dt = c_CFL*dt_CFL;
	
	// Restore vectors
	ierr = DMDAVecRestoreArrayDOF(db,localV,&v);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(da,localU,&u);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayDOF(da,F,&f);CHKERRQ(ierr);
	
	ierr = DMRestoreLocalVector(db,&localV);CHKERRQ(ierr);
	ierr = DMRestoreLocalVector(da,&localU);CHKERRQ(ierr);
	
	ierr = DMRestoreGlobalVector(db,&V); CHKERRQ(ierr);
	PetscFunctionReturn(0); 
} 
