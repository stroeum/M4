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
  PetscBool      Btype = user->BfieldType;
//PetscInt       mx = user->mx, my = user->my, mz = user->mz;
  PetscReal      *x=user->x, *y=user->y, *z=user->z;
  PetscReal      L=user->L;
  PetscReal      Xmin=user->outXmin, Ymin=user->outYmin, Zmin=user->outZmin;
  PetscReal      inZmax=user->inZmax;
//PetscReal      /*me = user->me,*/ mi[3] = {user->mi[O2p], user->mi[CO2p], user->mi[Op]};
  PetscReal      n0=user->n0, v0=user->v0, p0=user->p0, B0=user->B0;
  PetscReal      Bo=user->B[0], a=user->B[1], b=user->B[2];
//PetscReal      un[3]={user->un[0],user->un[1],user->un[2]}; // neutral wind

  PetscInt       i,j,k,l,xs,ys,zs,xm,ym,zm;
  Vec            localU;
  PetscReal      ****u;
  PetscReal      X,Y,Z;
  PetscReal      Te,Ti;
  PetscReal      nio[3],neo,vio[3],pio[3],peo;
  PetscInt       rank;
  PetscViewer    fViewer;
  char           fName[PETSC_MAX_PATH_LEN];
  PetscBool      flag;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"Form Initial Conditions: ...\n");
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  if (user->isInputFile) {
    // Read binary file containing the solution //
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
    // Map U on local U //
    ierr = DMGetLocalVector(da,&localU);CHKERRQ(ierr);

    ierr = DMDAVecGetArrayDOF(da,localU,&u);CHKERRQ(ierr);

    // Get local grid boundaries //
    ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
    
    for (k=zs; k<zs+zm; k++) {
      Z = Zmin + z[k]*L;
      for (j=ys; j<ys+ym; j++) {
        Y = Ymin + y[j]*L;
        for (i=xs; i<xs+xm; i++) {
          X = Xmin + x[i]*L;

          Te = Profile(8,Z);
          Ti = Profile(9,Z);
          assert(Te>0);
          assert(Ti>0);
          neo = Profile(7,Z);
          //if(!(neo>0)) neo = eps.n*n0;
          assert(neo>0);

          nio[O2p]  = Partition(O2p ,Z)*neo;
          nio[CO2p] = Partition(CO2p,Z)*neo;
          nio[Op]   = Partition(Op  ,Z)*neo;
          for (l=0; l<3; l++) {
            if (!(nio[l] > 0)) nio[l] = eps.n*n0;
            vio[l] = 100.0;
            //vio[l] = PetscSqrtScalar(3*kB*Ti/mi[l]);
            pio[l] = nio[l]*kB*Ti;
            //pio[l] = Profile(13+l,Z);
            assert(pio[l]>0);
          }
          //veo = PetscSqrtScalar(3*kB*Te/me)/v0;
          peo = neo*kB*Te;
          //peo = Profile(12,Z);
          assert(peo>0);
          if(Btype==0) {
            u[k][j][i][d.B[0]] = Profile(0,Z)/B0;
            u[k][j][i][d.B[1]] = Profile(1,Z)/B0;
            u[k][j][i][d.B[2]] = Profile(2,Z)/B0;
          } else if (Btype==1) {
            u[k][j][i][d.B[0]] = 0.0/B0;
            u[k][j][i][d.B[1]] = 0.0/B0;
            u[k][j][i][d.B[2]] = Bo /B0;
          } else if (Btype==2) {
            u[k][j][i][d.B[0]] = 0.0/B0;
            u[k][j][i][d.B[1]] = Bo /B0;
            u[k][j][i][d.B[2]] = 0.0/B0;
          } else if (Btype==3) {
            u[k][j][i][d.B[0]] = V_Dipole(Bo,a,X,Y,Z,0)/B0;
            u[k][j][i][d.B[1]] = V_Dipole(Bo,a,X,Y,Z,1)/B0;
            u[k][j][i][d.B[2]] = V_Dipole(Bo,a,X,Y,Z,2)/B0;
          } else if (Btype==4) {
            u[k][j][i][d.B[0]] = H_Dipole(Bo,a,X,Y,Z,0)/B0;
            u[k][j][i][d.B[1]] = H_Dipole(Bo,a,X,Y,Z,1)/B0;
            u[k][j][i][d.B[2]] = H_Dipole(Bo,a,X,Y,Z,2)/B0;
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
          }
          for (l=0; l<3; l++) {
            u[k][j][i][d.ni[l]] = nio[l]/n0;
            if (Btype==2) {
              u[k][j][i][d.vi[l][0]] = 0.0   /v0;
              u[k][j][i][d.vi[l][1]] = 0.0   /v0;
              u[k][j][i][d.vi[l][2]] = vio[l]/v0;
            } else {
              u[k][j][i][d.vi[l][0]] = 0.0   /v0;
              u[k][j][i][d.vi[l][1]] = vio[l]/v0;
              u[k][j][i][d.vi[l][2]] = 0.0   /v0;
            }
            //vmax  = PetscMax(vmax,vio[l]/v0);
            u[k][j][i][d.pi[l]] = pio[l]/p0;
          }
          //vmax = PetscMax(vmax,veo/v0);
          u[k][j][i][d.pe] = peo/p0;
        }
      }
    }

    ierr = DMDAVecRestoreArrayDOF(da,localU,&u);CHKERRQ(ierr);
    ierr = DMLocalToGlobalBegin(da,localU,INSERT_VALUES,U);CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(da,localU,INSERT_VALUES,U);CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(da,&localU);CHKERRQ(ierr);
  }

  PetscPrintf(PETSC_COMM_WORLD,"ARE WE THERE YET?\n");
  // Compute initial BC for u //
  ierr = FormBCu(U,user);CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD,"Form Initial Conditions: DONE\n");
  PetscFunctionReturn(0); 
} 

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "FormBCu"
/* 
   FormBoundaryConditions - Evaluates nonlinear function @ BC (outside the domain), F(u).
   * U: Up, D: Down, N: North, S: South, E: East, W: West
   * 8 corners: NWD, SWD, NED, SED, NWU, SWU, NEU, SEU
   * 12 edges: WD, ED, WU, EU, ND, SD, NU, SU, NE, NW, SE, SW
   * 6 faces: N, S, W, E, D, U
   Input Parameters:
.  u - input pointer

   Output Parameter:
.  u - output pointer
.  ctx - optional user-defined context, as set by SNESSetFunction()
*/
PetscErrorCode FormBCu(Vec U, void*ctx)
{
  PetscErrorCode ierr;
  AppCtx         *user = (AppCtx*)ctx;
  dvi            d = user->d;
  DM             da = user->da;
  PetscInt       mx = user->mx, my = user->my, mz = user->mz;
  PetscReal      *x=user->x,*y=user->y,*z=user->z;
  PetscInt       bcType = user->bcType;

  PetscReal      ****u;
  Vec            localU;
//PetscReal      vi[3][3],ve[3]; // velocity of i,e
  PetscInt       i,j,k,l,m,xs,ys,zs,xm,ym,zm;
  PetscInt       rank;
  coeff2         dh;
  coeff3         D1,D2;
  stencil        id;
  
  PetscFunctionBegin;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  // Get local grid boundaries //
  ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

  // Map U on localU //
  ierr = DMGetLocalVector(da,&localU);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(da,U,INSERT_VALUES,localU);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da,U,INSERT_VALUES,localU);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(da,localU,&u);CHKERRQ(ierr);


  for (k=zs; k<zs+zm; k++) {
    for (j=ys; j<ys+ym; j++) {
      if(xs==0) {                     // S-boundary //
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
      }
      if(xs+xm==mx) {                 // N-boundary //
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
      }
    }
  }
  for (i=xs; i<xs+xm; i++) {
    for (k=zs; k<zs+zm; k++) {
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
      }
    }
  }
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      if(zs==0) {                     // D-boundary //
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
      }
      if(zs+zm==mz) {                 // U-boundary //
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
      }
    }
  }
  /*
  if (rank == 0) { 
    PetscPrintf(PETSC_COMM_WORLD,"step = %d\n",user->cnt);
    for (k=-1;k<=1;k++) for (j=-1;j<=1;j++) for (i=-1;i<=1;i++) 
      if (! ( (i==-1 && j==-1) || (i==-1 && k==-1) || (j==-1 && k==-1) ) ) {
        PetscPrintf(PETSC_COMM_SELF,"@[%2d][%2d][%2d] nO2+=%2.5e; nCO2+=%2.5e; nO+=%2.5e; vO2+=[%2.5e %2.5e %2.5e]; vCO2+=[%2.5e %2.5e %2.5e]; vO+=[%2.5e %2.5e %2.5e]; pO2+=%2.5e; pCO2+=%2.5e; pO+=%2.5e; pe=%2.5e; B=[%2.5e %2.5e %2.5e]\n",i,j,k,u[k][j][i][0]*user->n0,u[k][j][i][1]*user->n0,u[k][j][i][2]*user->n0,u[k][j][i][3]*user->v0,u[k][j][i][4]*user->v0,u[k][j][i][5]*user->v0,u[k][j][i][6]*user->v0,u[k][j][i][7]*user->v0,u[k][j][i][8]*user->v0,u[k][j][i][9]*user->v0,u[k][j][i][10]*user->v0,u[k][j][i][11]*user->v0,u[k][j][i][12]*user->p0,u[k][j][i][13]*user->p0,u[k][j][i][14]*user->p0,u[k][j][i][15]*user->p0,u[k][j][i][16]*user->B0,u[k][j][i][17]*user->B0,u[k][j][i][18]*user->B0);
    }
  }
  */
  ierr = DMDAVecRestoreArrayDOF(da,localU,&u);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(da,localU,INSERT_VALUES,U);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(da,localU,INSERT_VALUES,U);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(da,&localU);CHKERRQ(ierr);
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

PetscErrorCode CheckValues(Vec U, Vec V,void *ctx)
{
  PetscErrorCode ierr;
  AppCtx         *user = (AppCtx*)ctx;
  DM             da = (DM)user->da;
  DM             db = (DM)user->db;
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
  PetscReal      ****u,****v;
  Vec            localU,localV;

  PetscFunctionBegin;
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  // Get local grid boundaries //
  ierr = DMDAGetGhostCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

  // Map U,V on localU, localV //
  ierr = DMGetLocalVector(da,&localU);CHKERRQ(ierr);
  ierr = DMGetLocalVector(db,&localV);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(da,U,INSERT_VALUES,localU);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(db,V,INSERT_VALUES,localV);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da,U,INSERT_VALUES,localU);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(db,V,INSERT_VALUES,localV);CHKERRQ(ierr);

  ierr = DMDAVecGetArrayDOF(da,localU,&u);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(db,localV,&v);CHKERRQ(ierr);
  
  // Check values //
  for (k=zs; k<zs+zm; k++) {
    Z = Zmin + z[k]*L; 
    for (j=ys; j<ys+ym; j++) {
      for (i=xs; i<xs+xm; i++) {
        if (limiters) {
          if(!(u[k][j][i][d.pe]>=eps.p)) {
            u[k][j][i][d.pe]=eps.p;
            //PetscPrintf(PETSC_COMM_SELF, "pe=%2.3e @ [%d,%d,%d] (Z=%f km)\n",u[k][j][i][d.pe]*p0, i,j,k, Z*1e-3);
          }
          assert(u[k][j][i][d.pe]>=eps.p); // check pe>=0
          for (l=0; l<3; l++) {
            if(!(u[k][j][i][d.ni[l]]>=eps.n)) {
              u[k][j][i][d.ni[l]]=eps.n;
              //PetscPrintf(PETSC_COMM_SELF, "ni[%d]=%2.3e @ [%d,%d,%d] (Z=%f km)\n",l,u[k][j][i][d.ni[l]]*n0, i,j,k, Z*1e-3);
            }
            assert(u[k][j][i][d.ni[l]]>=eps.n); // check ni>0
            if(!(u[k][j][i][d.pi[l]]>=eps.p)) {
              u[k][j][i][d.pi[l]]=eps.p;
              //PetscPrintf(PETSC_COMM_SELF, "pi[%d]=%2.3e @ [%d,%d,%d] (Z=%f km)\n",l,u[k][j][i][d.pi[l]]*p0, i,j,k, Z*1e-3);
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

  // Clean temporary variables //
  ierr = DMDAVecRestoreArrayDOF(db,localV,&v);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(da,localU,&u);CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(db,localV,INSERT_VALUES,V);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(da,localU,INSERT_VALUES,U);CHKERRQ(ierr);
  
  ierr = DMLocalToGlobalEnd(da,localU,INSERT_VALUES,U);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(da,&localU);CHKERRQ(ierr);

  ierr = DMLocalToGlobalEnd(db,localV,INSERT_VALUES,V);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(db,&localV);CHKERRQ(ierr);

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
PetscErrorCode FormIntermediateFunction(Vec U, Vec V,void *ctx)
{
  PetscErrorCode ierr;
  AppCtx         *user = (AppCtx*)ctx;
//tolerance      eps = user->eps;
  DM             da = (DM)user->da;
  DM             db = (DM)user->db;
  svi            s = user->s;
  
  PetscInt       mx = user->mx, my = user->my, mz = user->mz;
  PetscReal      n0 = user->n0, v0 = user->v0, p0 = user->p0;
  PetscReal      un[3] = {user->un[0]/v0, user->un[1]/v0, user->un[2]/v0};       // neutral wind
  PetscInt       rank;
  PetscInt       i,j,k,l,m,xs,ys,zs,xm,ym,zm;
  PetscReal      ****u,****v;
  Vec            localU,localV;
  PetscReal      vi[3][3],ve[3]; // velocity of i,e
  PetscReal      E[3];  // E-field
  PetscReal      J[3];  // total conduction current
  PetscReal      nu_en[2] = {0.0, 0.0};  // electron-neutral collisions

  dvi            d = user->d;
  PetscReal      L = user->L;
  PetscReal      Zmin = user->outZmin, Z;
  PetscReal      *x=user->x,*y=user->y,*z=user->z;
  PetscReal      tau = user->tau; 
  PetscReal      ni[3],ne,nn[2]; // density of i,e
  PetscReal      pe;
  PetscReal      Te;
  PetscReal      B[3];
  diff           dB[3],dpe;
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

  // Get local grid boundaries //
  ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

//printf("@%d, xs=%d\txs+xm=%d\tys=%d\tys+ym=%d\tzs=%d\tzs+zm=%d\n",rank,xs,xs+xm,ys,ys+ym,zs,zs+zm);
//exit(1);

  // Map U on local U //
  ierr = DMGetLocalVector(da,&localU);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(da,U,INSERT_VALUES,localU);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da,U,INSERT_VALUES,localU);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(da,localU,&u);CHKERRQ(ierr);

  // Create localV //
  ierr = DMGetLocalVector(db,&localV);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(db,localV,&v);CHKERRQ(ierr);
  /*
  if (rank == 0) { 
    printf("step == %d\n",user->cnt);
    i=0;
    j=0;
    k=0;
    PetscPrintf(PETSC_COMM_SELF,"@[%2d][%2d][%2d] B=[%2.5e %2.5e %2.5e];\n",i,j,k,u[k][j][i][d.B[0]]*user->B0,u[k-1][j][i][d.B[0]]*user->B0, u[k+1][j][i][d.B[0]]*user->B0);
  }
  */
 
  for (k=zs; k<zs+zm; k++) {
    Z = Zmin + z[k]*L;
    for (j=ys; j<ys+ym; j++) {
      for (i=xs; i<xs+xm; i++) {
        if (bcType == 0 || bcType == 1 || bcType == 2 || bcType == 10) {
          // Defines upper and lower indices for the differenciation //
          id.x[0] = i; id.x[1] = i-1; id.x[2] = i+1; 
          id.y[0] = j; id.y[1] = j-1; id.y[2] = j+1;
          id.z[0] = k; id.z[1] = k-1; id.z[2] = k+1;

          // Defines helpful coefficients for the differenciation //
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
          // Defines upper and lower indices for the differenciation //
          if (i==0)         {id.x[0] = i; id.x[1] = i+1; id.x[2] = i+2;}
          else if (i==mx-1) {id.x[0] = i; id.x[1] = i-2; id.x[2] = i-1;}
          else              {id.x[0] = i; id.x[1] = i-1; id.x[2] = i+1;}
          if (j==0)         {id.y[0] = j; id.y[1] = j+1; id.y[2] = j+2;}
          else if (j==my-1) {id.y[0] = j; id.y[1] = j-2; id.y[2] = j-1;}
          else              {id.y[0] = j; id.y[1] = j-1; id.y[2] = j+1;}
          if (k==0)         {id.z[0] = k; id.z[1] = k+1; id.z[2] = k+2;}
          else if (k==mz-1) {id.z[0] = k; id.z[1] = k-2; id.z[2] = k-1;}
          else              {id.z[0] = k; id.z[1] = k-1; id.z[2] = k+1;}

          // Defines helpful coefficients for the differenciation //
          dh.x[0] = x[id.x[1]] - x[id.x[0]]; dh.x[1] = x[id.x[2]] - x[id.x[0]];
          dh.y[0] = y[id.y[1]] - y[id.y[0]]; dh.y[1] = y[id.y[2]] - y[id.y[0]];
          dh.z[0] = z[id.z[1]] - z[id.z[0]]; dh.z[1] = z[id.z[2]] - z[id.z[0]];
        }

        // Calculate v on all interior points applies to all BC cases for v //
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

        // Calculate ne, ve, vi, J, E //
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

        Te          = pe*p0/(ne*n0*kB);

        nn[CO2]     = Profile(5,Z);
        nn[O]       = Profile(4,Z);

        // Intermediate calculations //
        dpe.dx = D1.x[0]*u[k][j][id.x[0]][d.pe] + D1.x[1]*u[k][j][id.x[1]][d.pe] + D1.x[2]*u[k][j][id.x[2]][d.pe]; //dpe.dx 
        dpe.dy = D1.y[0]*u[k][id.y[0]][i][d.pe] + D1.y[1]*u[k][id.y[1]][i][d.pe] + D1.y[2]*u[k][id.y[2]][i][d.pe]; //dpe.dy 
        dpe.dz = D1.z[0]*u[id.z[0]][j][i][d.pe] + D1.z[1]*u[id.z[1]][j][i][d.pe] + D1.z[2]*u[id.z[2]][j][i][d.pe]; //dpe.dz 

        for (m=0; m<3; m++) {
          //dB[m].dz = (u[idx.z[2]][j][i][d.B[m]]-u[idx.z[0]][j][i][d.B[m]])/(z[idx.z[2]]-z[idx.z[0]]); //dBy.dz 
          dB[m].dx = D1.x[0]*u[k][j][id.x[0]][d.B[m]] + D1.x[1]*u[k][j][id.x[1]][d.B[m]] + D1.x[2]*u[k][j][id.x[2]][d.B[m]]; //dB[m].dx 
          dB[m].dy = D1.y[0]*u[k][id.y[0]][i][d.B[m]] + D1.y[1]*u[k][id.y[1]][i][d.B[m]] + D1.y[2]*u[k][id.y[2]][i][d.B[m]]; //dB[m].dy 
          dB[m].dz = D1.z[0]*u[id.z[0]][j][i][d.B[m]] + D1.z[1]*u[id.z[1]][j][i][d.B[m]] + D1.z[2]*u[id.z[2]][j][i][d.B[m]]; //dB[m].dz 
        }

        nu_en[CO2] = Ven(CO2, nn[CO2], Te);
        nu_en[O]   = Ven(O  , nn[O]  , Te);

        // Current //
        J[0] = dB[2].dy - dB[1].dz; 
        J[1] = dB[0].dz - dB[2].dx;
        J[2] = dB[1].dx - dB[0].dy;
       
        // Electron velocity //
        ve[0] = ( (ni[O2p]*vi[O2p][0] + ni[CO2p]*vi[CO2p][0] + ni[Op]*vi[Op][0]) - J[0] )/ ne;
        ve[1] = ( (ni[O2p]*vi[O2p][1] + ni[CO2p]*vi[CO2p][1] + ni[Op]*vi[Op][1]) - J[1] )/ ne;
        ve[2] = ( (ni[O2p]*vi[O2p][2] + ni[CO2p]*vi[CO2p][2] + ni[Op]*vi[Op][2]) - J[2] )/ ne;
        //if (i==0 && j==0 && k==0) printf("dz=%2.1e km, nu_en=%2.5e ve=[%2.3e %2.3e %2.3e] un=[%2.3e %2.3e %2.3e]\n", dh.z[0],nu_en[CO2]+nu_en[O], ve[0]*v0,ve[1]*v0,ve[2]*v0,un[0]*v0,un[1]*v0,un[2]*v0);
        
        // E-field //
        E[0] = -CrossP(ve,B,0) - dpe.dx/ne + (nu_en[CO2] + nu_en[O]) * tau * (un[0] - ve[0]);
        E[1] = -CrossP(ve,B,1) - dpe.dy/ne + (nu_en[CO2] + nu_en[O]) * tau * (un[1] - ve[1]);
        E[2] = -CrossP(ve,B,2) - dpe.dz/ne + (nu_en[CO2] + nu_en[O]) * tau * (un[2] - ve[2]);

        // Store data //
        for (m=0; m<3; m++) { 
          v[k][j][i][s.ve[m]] = ve[m];
          v[k][j][i][s.J[m]]  = J[m];
          v[k][j][i][s.E[m]]  = E[m];
        }
      }
    }
  }
  // Calculate BC points for V //
  for (j=ys; j<ys+ym; j++) {
    for (i=xs; i<xs+xm; i++) {
      if(zs==0) {                     // D-boundary //
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
      if(zs+zm==mz) {                 // U-boundary //
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
      if(xs==0) {                     // S-boundary //
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
      if(xs+xm==mx) {                 // N-boundary //
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
  /*
  if (rank == 7) { 
    for (k=mz-2;k<=mz;k++) for (j=my-2;j<=my;j++) for (i=mx-2;i<=mx;i++) 
      if (! ( (i==mx && j==my) || (i==mx && k==mz) || (j==my && k==mz) ) ) {
        PetscPrintf(PETSC_COMM_SELF,"@[%2d][%2d][%2d] ve=[%2.5e %2.5e %2.5e]; J=[%2.5e %2.5e %2.5e]; E=[%2.5e %2.5e %2.5e]\n",i,j,k,v[k][j][i][0]*user->v0,v[k][j][i][1]*user->v0,v[k][j][i][2]*user->v0,v[k][j][i][3]*qe*user->n0*user->v0,v[k][j][i][4]*qe*user->n0*user->v0,v[k][j][i][5]*qe*user->n0*user->v0,v[k][j][i][6]*user->v0*user->B0,v[k][j][i][7]*user->v0*user->B0,v[k][j][i][8]*user->v0*user->B0);
      }
    }
    */

  ierr = DMDAVecRestoreArrayDOF(db,localV,&v);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(da,localU,&u);CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(db,localV,INSERT_VALUES,V);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(da,localU,INSERT_VALUES,U);CHKERRQ(ierr);
  
  ierr = DMLocalToGlobalEnd(da,localU,INSERT_VALUES,U);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(da,&localU);CHKERRQ(ierr);

  ierr = DMLocalToGlobalEnd(db,localV,INSERT_VALUES,V);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(db,&localV);CHKERRQ(ierr);
 
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
  PetscReal      Zmin = user->outZmin;
  PetscReal      n0 = user->n0, v0 = user->v0, p0 = user->p0;
  DM             da = (DM)user->da;
  DM             db = (DM)user->db;
  PetscInt       mx = user->mx, my = user->my, mz = user->mz;
  PetscReal      L=user->L;
  PetscReal      *x=user->x,*y=user->y,*z=user->z;
  PetscInt       bcType = user->bcType;
  dvi            d = user->d;
  svi            s = user->s;
  PetscReal      rM = user->rM;
  PetscReal      tau = user->tau; 
  PetscReal      me = user->me;                                                  // me 
  PetscReal      mi[3] = {user->mi[0], user->mi[1], user->mi[2]};                // mi[O2+,CO2+,O+]
  PetscReal      un[3] = {user->un[0]/v0, user->un[1]/v0, user->un[2]/v0};       // neutral wind
  PetscReal      gama[4] = {user->gama[O2p],user->gama[CO2p],user->gama[Op],user->gama[e]};
  PetscReal      vp      = 0.0;                                                  // Alfven speed
  PetscReal      vCFL[3] = {0.0, 0.0, 0.0};                                      // x-, y-, and z-CFL criterion v
  PetscReal      tA  =0.0, vA  =0.0;                                             // Alfven speed
  PetscReal      ts_n=0.0, vs_n=0.0;                                             // speed of sound
  PetscReal      ts_c=0.0, vs_c=0.0;                                             // speed of sound
  PetscReal      rhoc=0.0;                                                       // kg/m^3, weighted mass density of charged species (i,e)
  PetscErrorCode ierr;
  PetscReal      dt,dtCFL,h[3];
  PetscReal      nu[4][2];                                                       // ion-neutral collision frequency
  PetscReal      E[3],B[3];                                                      // E-field, B-field
  PetscReal      ne,ni[3],nn[2];                                                 // [O2+], [CO2+], [O+], ne in m^-3
  PetscReal      Ti[3],Te,Tn;                                                    // Temperature of e, O2+, CO2+, O+, CO2, 0 in K
  PetscReal      pe,pi[3];                                                       // pressure of O2+, CO2+, O+, e in J/m^3
  PetscReal      ****u,****v,****f;
  Vec            V;
  Vec            localU,localV;
  PetscInt       rank;
  PetscInt       i,j,k,l,m,n,xs,ys,zs,xm,ym,zm;
  PetscReal      Z;
  PetscReal      ve[3],vi[3][3];                                                 // electron and ion velocities 
  PetscReal      ve_max[3],vi_max[3][3];
  diff           dE[3],dpe,dve[3],dni[3],dpi[3],dvi[3][3];                       // Other derivatives 
  stencil        id;
  coeff3         D1;
  coeff2         dh;
  fdv            ff,ff_max,vf,vf_max;            
  PetscBool      flag = PETSC_FALSE;
  
  PetscFunctionBegin;
  id.x[0] = 0; id.y[0] = 0; id.z[0] = 0;
  id.x[1] = 0; id.y[1] = 0; id.z[1] = 0;
  id.x[2] = 0; id.y[2] = 0; id.z[2] = 0;

  dh.x[0] = 0; dh.y[0] = 0; dh.z[0] = 0;
  dh.x[1] = 0; dh.y[1] = 0; dh.z[1] = 0;
  
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  ierr = TSGetTimeStep(ts,&dt);CHKERRQ(ierr);
  ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  
  // Form boundary conditions //
  //PetscPrintf(PETSC_COMM_WORLD,"Form Boundary Conditions: ...\n");
  ierr = FormBCu(U,user);CHKERRQ(ierr);
  //PetscPrintf(PETSC_COMM_WORLD,"Form Boundary Conditions: DONE\n");

  // Intermediate calculations //
  ierr = DMGetGlobalVector(db,&V);CHKERRQ(ierr);
  //PetscPrintf(PETSC_COMM_WORLD,"Form Intermediate Function: ...\n");
  ierr = FormIntermediateFunction(U,V,user);CHKERRQ(ierr);
  //PetscPrintf(PETSC_COMM_WORLD,"Form Intermediate Function: DONE\n");

  // Check values sanity //
  //PetscPrintf(PETSC_COMM_WORLD,"Check Values Sanity: ...\n");
  ierr = CheckValues(U,V,user);CHKERRQ(ierr);
  //PetscPrintf(PETSC_COMM_WORLD,"Check Values Sanity: DONE\n");

  // Map U,V on localU,localV //
  ierr = DMGetLocalVector(da,&localU);CHKERRQ(ierr);
  ierr = DMGetLocalVector(db,&localV);CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(da,U,INSERT_VALUES,localU);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(db,V,INSERT_VALUES,localV);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(db,V,INSERT_VALUES,localV);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da,U,INSERT_VALUES,localU);CHKERRQ(ierr);

  ierr = DMDAVecGetArrayDOF(da,localU,&u);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(db,localV,&v);CHKERRQ(ierr);
  ierr = DMDAVecGetArrayDOF(da,F,&f);CHKERRQ(ierr);

  // Initialize max speed //

  for (m=0; m<3; m++) {
    ve_max[m] = 0.0;
    for (l=0; l<3; l++) {
      vi_max[l][m]    = 0.0;
    }
  }
  for (m=0; m<3; m++) {
    vf_max.g[m] = 0.0;
    ff_max.g[m] = 0.0;
    for (l=0; l<4; l++) {
      vf_max.dp[l][m]  = 0.0;
      vf_max.E[l][m]   = 0.0;
      vf_max.B[l][m]   = 0.0;
      vf_max.c[l][m]   = 0.0;
      vf_max.adv[l][m] = 0.0;
      ff_max.dp[l][m]  = 0.0;
      ff_max.E[l][m]   = 0.0;
      ff_max.B[l][m]   = 0.0;
      ff_max.c[l][m]   = 0.0;
      ff_max.adv[l][m] = 0.0;
    }
  }

  // Form F(u,v) = Udot //
  
  // Compute function over the locally owned part of the grid //
  dtCFL = 1e300;
  flag = PETSC_FALSE;
  for (k=zs; k<zs+zm; k++) {
    Z = Zmin + z[k]*L;
    for (j=ys; j<ys+ym; j++) {
      for (i=xs; i<xs+xm; i++) {
        if (bcType == 0 || bcType == 1 || bcType == 2 || bcType == 10) {
          // Defines upper and lower indices for the differenciation //
          id.x[0] = i; id.x[1] = i-1; id.x[2] = i+1; 
          id.y[0] = j; id.y[1] = j-1; id.y[2] = j+1;
          id.z[0] = k; id.z[1] = k-1; id.z[2] = k+1;

          // Defines helpful coefficients for the differenciation //
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
          // Defines upper and lower indices for the differenciation //
          if (i==0)         {id.x[0] = i; id.x[1] = i+1; id.x[2] = i+2;}
          else if (i==mx-1) {id.x[0] = i; id.x[1] = i-2; id.x[2] = i-1;}
          else              {id.x[0] = i; id.x[1] = i-1; id.x[2] = i+1;}
          if (j==0)         {id.y[0] = j; id.y[1] = j+1; id.y[2] = j+2;}
          else if (j==my-1) {id.y[0] = j; id.y[1] = j-2; id.y[2] = j-1;}
          else              {id.y[0] = j; id.y[1] = j-1; id.y[2] = j+1;}
          if (k==0)         {id.z[0] = k; id.z[1] = k+1; id.z[2] = k+2;}
          else if (k==mz-1) {id.z[0] = k; id.z[1] = k-2; id.z[2] = k-1;}
          else              {id.z[0] = k; id.z[1] = k-1; id.z[2] = k+1;}

          // Defines helpful coefficients for the differenciation //
          dh.x[0] = x[id.x[1]] - x[id.x[0]]; dh.x[1] = x[id.x[2]] - x[id.x[0]];
          dh.y[0] = y[id.y[1]] - y[id.y[0]]; dh.y[1] = y[id.y[2]] - y[id.y[0]];
          dh.z[0] = z[id.z[1]] - z[id.z[0]]; dh.z[1] = z[id.z[2]] - z[id.z[0]];
        }

        // Coefficients for differenciation //
        D1.x[0]  = - 1.0/dh.x[0] - 1.0/dh.x[1]; 
        D1.x[1]  =   1.0/dh.x[0] - 1.0/(dh.x[0]-dh.x[1]); 
        D1.x[2]  =   1.0/dh.x[1] + 1.0/(dh.x[0]-dh.x[1]);
        D1.y[0]  = - 1.0/dh.y[0] - 1.0/dh.y[1]; 
        D1.y[1]  =   1.0/dh.y[0] - 1.0/(dh.y[0]-dh.y[1]); 
        D1.y[2]  =   1.0/dh.y[1] + 1.0/(dh.y[0]-dh.y[1]);
        D1.z[0]  = - 1.0/dh.z[0] - 1.0/dh.z[1]; 
        D1.z[1]  =   1.0/dh.z[0] - 1.0/(dh.z[0]-dh.z[1]); 
        D1.z[2]  =   1.0/dh.z[1] + 1.0/(dh.z[0]-dh.z[1]);

        // For readability, we distribute the components of u onto the variable n, v, p, B //
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

        nn[CO2] = Profile(5,Z); nn[O] = Profile(4,Z);
        Tn = Profile(10,Z);
        
        for (l=0; l<3; l++) {
          Ti[l] = pi[l]*p0/(ni[l]*n0*kB);
          for (n=0; n<2; n++) {
            nu[l][n] = Vin(l,n,nn[n],Ti[l],Tn);
          }
        }
        Te = pe*p0/(ne*n0*kB);
        for (n=0; n<2; n++) 
          nu[e][n] = Ven(n,nn[n],Te);

        //if (i==mx-1 && j==my-1 && k==mz-1) printf("nu.O2pn=%+14.7e; nu.CO2pn=%+14.7e; nu.Opn=%+14.7e; nu.en=%+14.7e;\n",nu[O2p][CO2]+nu[O2p][O],nu[CO2p][CO2]+nu[CO2p][O],nu[Op][CO2]+nu[Op][O],nu[e][CO2]+nu[e][O]);

        // Intermediate calculations //
        for (m=0; m<3; m++) {
          dE[m].dx = D1.x[0]*v[k][j][id.x[0]][s.E[m]] + D1.x[1]*v[k][j][id.x[1]][s.E[m]] + D1.x[2]*v[k][j][id.x[2]][s.E[m]]; //dE[m].dx 
          dE[m].dy = D1.y[0]*v[k][id.y[0]][i][s.E[m]] + D1.y[1]*v[k][id.y[1]][i][s.E[m]] + D1.y[2]*v[k][id.y[2]][i][s.E[m]]; //dE[m].dy 
          dE[m].dz = D1.z[0]*v[id.z[0]][j][i][s.E[m]] + D1.z[1]*v[id.z[1]][j][i][s.E[m]] + D1.z[2]*v[id.z[2]][j][i][s.E[m]]; //dE[m].dz 

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

        // Equations //
        rhoc = ni[O2p]*mi[O2p]/me + ni[CO2p]*mi[CO2p]/me + ni[Op]*mi[Op]/me + ne;
        vs_n = Profile(16,Z)/v0;
        vs_c = PetscSqrtScalar((gama[O2p]*pi[O2p]+gama[CO2p]*pi[CO2p]+gama[Op]*pi[Op] + gama[e]*pe)/rhoc);
        vA   = PetscSqrtScalar(Norm2(B)/rhoc);

        // CFL // 
        vCFL[0]=0.0; vCFL[1]=0.0; vCFL[2]=0.0;

        // gravity force //
        ff.g[0] =   0.0;
        ff.g[1] =   0.0;
        ff.g[2] = - 1.0/((1+Z/rM)*(1+Z/rM));

        for (l=0; l<3; l++) {
          // grap-pi force //
          ff.dp[l][0] = -me/mi[l]*dpi[l].dx/ni[l];
          ff.dp[l][1] = -me/mi[l]*dpi[l].dy/ni[l];
          ff.dp[l][2] = -me/mi[l]*dpi[l].dz/ni[l]; 
          // advection //
          ff.adv[l][0] = - (vi[l][0]*dvi[l][0].dx +vi[l][1]*dvi[l][0].dy +vi[l][2]*dvi[l][0].dz);
          ff.adv[l][1] = - (vi[l][0]*dvi[l][1].dx +vi[l][1]*dvi[l][1].dy +vi[l][2]*dvi[l][1].dz);
          ff.adv[l][2] = - (vi[l][0]*dvi[l][2].dx +vi[l][1]*dvi[l][2].dy +vi[l][2]*dvi[l][2].dz);
          for (m=0; m<3; m++) {
            // E-Lorentz force //
            ff.E[l][m] =  me/mi[l]*E[m];
            // B-Lorentz force //
            ff.B[l][m] =  me/mi[l]*CrossP(vi[l],B,m);
            // collisions //
            //ff.c[l][m] = tau*(nu[l][CO2]+nu[l][O])*(PetscAbsScalar(un[m])+PetscAbsScalar(vi[l][m]));
            ff.c[l][m] = tau*(nu[l][CO2]+nu[l][O])*(un[m]-vi[l][m]);
          }
        }

        // grad-pe force //
        ff.dp[e][0] = -dpe.dx/ne;
        ff.dp[e][1] = -dpe.dy/ne;
        ff.dp[e][2] = -dpe.dz/ne;  
        ff.adv[e][0] = - (ve[0]*dve[0].dx +ve[1]*dve[0].dy +ve[2]*dve[0].dz);
        ff.adv[e][1] = - (ve[0]*dve[1].dx +ve[1]*dve[1].dy +ve[2]*dve[1].dz);
        ff.adv[e][2] = - (ve[0]*dve[2].dx +ve[1]*dve[2].dy +ve[2]*dve[2].dz);

        for (m=0; m<3; m++) {
          // E-Lorentz force //
          ff.E[e][m] = -E[m];
          // B-Lorentz force //
          ff.B[e][m] =  -CrossP(ve,B,m);
          // collisions //
          //ff.c[e][m] = tau*(nu[e][CO2]+nu[e][O])*(PetscAbsScalar(un[m])+PetscAbsScalar(ve[m]));
          ff.c[e][m] = tau*(nu[e][CO2]+nu[e][O])*(un[m]-ve[m]);
        }

        for (m=0; m<3; m++) {
          vCFL[m] = MaxAbs(vCFL[m],vA  );
          vCFL[m] = MaxAbs(vCFL[m],vs_n);
          vCFL[m] = MaxAbs(vCFL[m],vs_c);   
          vCFL[m] = MaxAbs(vCFL[m],ve[m]); ve_max[m] = MaxAbs(ve_max[m],ve[m]);
          for (l=0; l<3; l++) {vCFL[m] = MaxAbs(vCFL[m],vi[l][m]); vi_max[l][m] = MaxAbs(vi_max[l][m],vi[l][m]);}
        }

        for (m=0; m<3; m++) {
          vf.g[m] = CrossP(ff.g,B,m)/Norm2(B); vf_max.g[m] = MaxAbs(vf_max.g[m],PetscAbsScalar(vf.g[m])); vCFL[m] = MaxAbs(vCFL[m],PetscAbsScalar(vf.g[m])); 
          ff_max.g[m]  = MaxAbs(ff_max.g[m] ,PetscAbsScalar(ff.g[m]));

          for (l=0; l<4; l++) {
            vf.adv[l][m] = CrossP(ff.adv[l],B,m)/Norm2(B); vf_max.adv[l][m] = MaxAbs(vf_max.adv[l][m],vf.adv[l][m]); vCFL[m] = MaxAbs(vCFL[m],vf.adv[l][m]); 
            vf.dp[l][m]  = CrossP(ff.dp[l] ,B,m)/Norm2(B); vf_max.dp[l][m]  = MaxAbs(vf_max.dp[l][m] ,vf.dp[l][m] ); vCFL[m] = MaxAbs(vCFL[m],vf.dp[l][m] ); 
            vf.E[l][m]   = CrossP(ff.E[l]  ,B,m)/Norm2(B); vf_max.E[l][m]   = MaxAbs(vf_max.E[l][m]  ,vf.E[l][m]  ); vCFL[m] = MaxAbs(vCFL[m],vf.E[l][m]  ); 
            vf.B[l][m]   = CrossP(ff.B[l]  ,B,m)/Norm2(B); vf_max.B[l][m]   = MaxAbs(vf_max.B[l][m]  ,vf.B[l][m]  ); vCFL[m] = MaxAbs(vCFL[m],vf.B[l][m]  ); 
            vf.c[l][m]   = CrossP(ff.c[l]  ,B,m)/Norm2(B); vf_max.c[l][m]   = MaxAbs(vf_max.c[l][m]  ,vf.c[l][m]  ); vCFL[m] = MaxAbs(vCFL[m],vf.c[l][m]  ); 
            ff_max.adv[l][m] = MaxAbs(ff_max.adv[l][m],ff.adv[l][m]);
            ff_max.dp[l][m]  = MaxAbs(ff_max.dp[l][m] ,ff.dp[l][m] );
            ff_max.E[l][m]   = MaxAbs(ff_max.E[l][m]  ,ff.E[l][m]  );
            ff_max.B[l][m]   = MaxAbs(ff_max.B[l][m]  ,ff.B[l][m]  );
            ff_max.c[l][m]   = MaxAbs(ff_max.c[l][m]  ,ff.c[l][m]  );
          }
        }
        h[0] = MinAbs(dh.x[0],dh.x[1]);
        h[1] = MinAbs(dh.y[0],dh.y[1]);
        h[2] = MinAbs(dh.z[0],dh.z[1]);
        if(i==0 && j==0) {
          tA   += h[2]/vA;
          ts_n += h[2]/vs_n;
          ts_c += h[2]/vs_c;
        } 

        if(dtCFL>1.0/(vCFL[0]/h[0]+vCFL[1]/h[1]+vCFL[2]/h[2]))
          dtCFL=1.0/(vCFL[0]/h[0]+vCFL[1]/h[1]+vCFL[2]/h[2]);
        
        vp = PetscSqrtScalar( MaxAbs( MaxAbs(Norm2(ve),Norm2(vi[O2p])) , MaxAbs(Norm2(vi[CO2p]),Norm2(vi[Op])) ) ); 
        /*
        if(vp > vA) PetscPrintf(PETSC_COMM_SELF, "SuperAlfvenic speed @ [%d,%d,%d], |vloc|=%2.3e, vA=%2.3e (Z=%f km)\n",i,j,k,vp*v0,vA*v0, Z*1e-3);
        assert(vp <= vA); // Check sub-Alfvenic speeds
        if(vp > vs_n) PetscPrintf(PETSC_COMM_SELF, "Supersonic (neutral) speed @ [%d,%d,%d], |vloc|=%2.3e, vs_n=%2.3e (Z=%f km)\n",i,j,k,vp*v0,vs_n*v0, Z*1e-3);
        assert(vp <= vs_n); // Check subsonic speeds
        if(vp > vs_c) PetscPrintf(PETSC_COMM_SELF, "Supersonic (charged) speed @ [%d,%d,%d], |vloc|=%2.3e, vs_c=%2.3e (Z=%f km)\n",i,j,k,vp*v0,vs_c*v0, Z*1e-3);
        assert(vp <= vs_c); // Check subsonic speeds
        */

        // Eq 0-2: Continuity equations for ions //
        for (l=0; l<3; l++)
          f[k][j][i][d.ni[l]] = -ni[l]*(dvi[l][0].dx + dvi[l][1].dy + dvi[l][2].dz) - (vi[l][0]*dni[l].dx + vi[l][1]*dni[l].dy + vi[l][2]*dni[l].dz);
      
        // Eq 3-11: Ion momenta // NEGLECT O and e-n collisions FOR NOW (07-29-11) //
        for (l=0; l<3; l++) {
          f[k][j][i][d.vi[l][0]] = - (vi[l][0]*dvi[l][0].dx +vi[l][1]*dvi[l][0].dy +vi[l][2]*dvi[l][0].dz) + me/mi[l]*(E[0] +CrossP(vi[l],B,0) -dpi[l].dx/ni[l])                           +tau*(nu[l][CO2]+nu[l][O])*(un[0]-vi[l][0]);
          f[k][j][i][d.vi[l][1]] = - (vi[l][0]*dvi[l][1].dx +vi[l][1]*dvi[l][1].dy +vi[l][2]*dvi[l][1].dz) + me/mi[l]*(E[1] +CrossP(vi[l],B,1) -dpi[l].dy/ni[l])                           +tau*(nu[l][CO2]+nu[l][O])*(un[1]-vi[l][1]);
          f[k][j][i][d.vi[l][2]] = - (vi[l][0]*dvi[l][2].dx +vi[l][1]*dvi[l][2].dy +vi[l][2]*dvi[l][2].dz) + me/mi[l]*(E[2] +CrossP(vi[l],B,2) -dpi[l].dz/ni[l]) - 1.0/((1+Z/rM)*(1+Z/rM)) +tau*(nu[l][CO2]+nu[l][O])*(un[2]-vi[l][2]); 
        }

        // Eq 12-14: Equations of state for ions //
        for (l=0; l<3; l++)
          f[k][j][i][d.pi[l]] = - gama[l] * pi[l] * (dvi[l][0].dx + dvi[l][1].dy + dvi[l][2].dz) - (vi[l][0]*dpi[l].dx + vi[l][1]*dpi[l].dy + vi[l][2]*dpi[l].dz);

        // Eq 15: Equation of state for electrons //
        f[k][j][i][d.pe]    = - gama[l] * pe    * (dve[0].dx    + dve[1].dy    + dve[2].dz   ) - (   ve[0]*dpe.dx    +    ve[1]*dpe.dy    +    ve[2]*dpe.dz   );
         
        // Eq 16-18: Faraday's law //
        f[k][j][i][d.B[0]] = - (dE[2].dy - dE[1].dz); 
        f[k][j][i][d.B[1]] = - (dE[0].dz - dE[2].dx);
        f[k][j][i][d.B[2]] = - (dE[1].dx - dE[0].dy);

        /*
        if(i==mx-1 && j==my-1 && k==mz-1) { 
          ierr = PetscPrintf(PETSC_COMM_SELF,"\tCONTRIBUTIONS: v.O2p     =[%+14.7e %+14.7e %+14.7e]; v.CO2p     =[%+14.7e %+14.7e %+14.7e]; v.Op     =[%+14.7e %+14.7e %+14.7e]; v.e     =[%+14.7e %+14.7e %+14.7e]\nff.g      =[%+14.7e %+14.7e %+14.7e]\nff.dp.O2p =[%+14.7e %+14.7e %+14.7e]; ff.dp.CO2p =[%+14.7e %+14.7e %+14.7e]; ff.dp.Op =[%+14.7e %+14.7e %+14.7e]; ff.dp.e =[%+14.7e %+14.7e %+14.7e]\nff.E.O2p  =[%+14.7e %+14.7e %+14.7e]; ff.E.CO2p  =[%+14.7e %+14.7e %+14.7e]; ff.E.Op  =[%+14.7e %+14.7e %+14.7e]; ff.E.e  =[%+14.7e %+14.7e %+14.7e]\nff.B.O2p  =[%+14.7e %+14.7e %+14.7e]; ff.B.CO2p  =[%+14.7e %+14.7e %+14.7e]; ff.B.Op  =[%+14.7e %+14.7e %+14.7e]; ff.B.e  =[%+14.7e %+14.7e %+14.7e]\nff.c.O2p  =[%+14.7e %+14.7e %+14.7e]; ff.c.CO2p  =[%+14.7e %+14.7e %+14.7e]; ff.c.Op  =[%+14.7e %+14.7e %+14.7e]; ff.c.e  =[%+14.7e %+14.7e %+14.7e]\nff.adv.O2p=[%+14.7e %+14.7e %+14.7e]; ff.adv.CO2p=[%+14.7e %+14.7e %+14.7e]; ff.adv.Op=[%+14.7e %+14.7e %+14.7e]; ff.adv.e=[%+14.7e %+14.7e %+14.7e]\nf.O2p     =[%+14.7e %+14.7e %+14.7e]; f.CO2p     =[%+14.7e %+14.7e %+14.7e]; f.Op     =[%+14.7e %+14.7e %+14.7e];\n", 
              vi[O2p][0]     *v0    ,vi[O2p][1]     *v0    ,vi[O2p][2]     *v0    ,
              vi[CO2p][0]    *v0    ,vi[CO2p][1]    *v0    ,vi[CO2p][2]    *v0    ,
              vi[Op][0]      *v0    ,vi[Op][1]      *v0    ,vi[Op][2]      *v0    ,
              ve[0]          *v0    ,ve[1]          *v0    ,ve[2]          *v0    ,
              ff.g[0]        *v0/tau,ff.g[1]        *v0/tau,ff.g[2]        *v0/tau,
              ff.dp[O2p][0]  *v0/tau,ff.dp[O2p][1]  *v0/tau,ff.dp[O2p][2]  *v0/tau,
              ff.dp[CO2p][0] *v0/tau,ff.dp[CO2p][1] *v0/tau,ff.dp[CO2p][2] *v0/tau,
              ff.dp[Op][0]   *v0/tau,ff.dp[Op][1]   *v0/tau,ff.dp[Op][2]   *v0/tau,
              ff.dp[e][0]    *v0/tau,ff.dp[e][1]    *v0/tau,ff.dp[e][2]    *v0/tau,
              ff.E[O2p][0]   *v0/tau,ff.E[O2p][1]   *v0/tau,ff.E[O2p][2]   *v0/tau,
              ff.E[CO2p][0]  *v0/tau,ff.E[CO2p][1]  *v0/tau,ff.E[CO2p][2]  *v0/tau,
              ff.E[Op][0]    *v0/tau,ff.E[Op][1]    *v0/tau,ff.E[Op][2]    *v0/tau,
              ff.E[e][0]     *v0/tau,ff.E[e][1]     *v0/tau,ff.E[e][2]     *v0/tau,
              ff.B[O2p][0]   *v0/tau,ff.B[O2p][1]   *v0/tau,ff.B[O2p][2]   *v0/tau,
              ff.B[CO2p][0]  *v0/tau,ff.B[CO2p][1]  *v0/tau,ff.B[CO2p][2]  *v0/tau,
              ff.B[Op][0]    *v0/tau,ff.B[Op][1]    *v0/tau,ff.B[Op][2]    *v0/tau,
              ff.B[e][0]     *v0/tau,ff.B[e][1]     *v0/tau,ff.B[e][2]     *v0/tau,
              ff.c[O2p][0]   *v0/tau,ff.c[O2p][1]   *v0/tau,ff.c[O2p][2]   *v0/tau,
              ff.c[CO2p][0]  *v0/tau,ff.c[CO2p][1]  *v0/tau,ff.c[CO2p][2]  *v0/tau,
              ff.c[Op][0]    *v0/tau,ff.c[Op][1]    *v0/tau,ff.c[Op][2]    *v0/tau,
              ff.c[e][0]     *v0/tau,ff.c[e][1]     *v0/tau,ff.c[e][2]     *v0/tau,
              ff.adv[O2p][0] *v0/tau,ff.adv[O2p][1] *v0/tau,ff.adv[O2p][2] *v0/tau,
              ff.adv[CO2p][0]*v0/tau,ff.adv[CO2p][1]*v0/tau,ff.adv[CO2p][2]*v0/tau,
              ff.adv[Op][0]  *v0/tau,ff.adv[Op][1]  *v0/tau,ff.adv[Op][2]  *v0/tau,
              ff.adv[e][0]   *v0/tau,ff.adv[e][1]   *v0/tau,ff.adv[e][2]   *v0/tau,
              f[k][j][i][d.vi[O2p][0]] *v0/tau,f[k][j][i][d.vi[O2p][1]] *v0/tau,f[k][j][i][d.vi[O2p][2]] *v0/tau,
              f[k][j][i][d.vi[CO2p][0]]*v0/tau,f[k][j][i][d.vi[CO2p][1]]*v0/tau,f[k][j][i][d.vi[CO2p][2]]*v0/tau,
              f[k][j][i][d.vi[Op][0]]  *v0/tau,f[k][j][i][d.vi[Op][1]]  *v0/tau,f[k][j][i][d.vi[Op][2]]  *v0/tau
                ); CHKERRQ(ierr);
        }
        */
        /*
        for (l=0; l<19; l++) 
          if (isnan(u[k][j][i][l])) 
            flag = PETSC_TRUE;
        if (flag) {
          PetscPrintf(PETSC_COMM_SELF,"@[%3d][%3d][%3d] nO2p=%+9.3e; nCO2p=%+9.3e; nOp=%+9.3e; vO2p=[%+9.3e %+9.3e %+9.3e]; vCO2p=[%+9.3e %+9.3e %+9.3e]; vOp=[%+9.3e %+9.3e %+9.3e]; pO2p=%+9.3e; pCO2p=%+9.3e; pOp=%+9.3e; pe=%+9.3e; B=[%+9.3e %+9.3e %+9.3e]\n",i,j,k,u[k][j][i][0]*user->n0,u[k][j][i][1]*user->n0,u[k][j][i][2]*user->n0,u[k][j][i][3]*user->v0,u[k][j][i][4]*user->v0,u[k][j][i][5]*user->v0,u[k][j][i][6]*user->v0,u[k][j][i][7]*user->v0,u[k][j][i][8]*user->v0,u[k][j][i][9]*user->v0,u[k][j][i][10]*user->v0,u[k][j][i][11]*user->v0,u[k][j][i][12]*user->p0,u[k][j][i][13]*user->p0,u[k][j][i][14]*user->p0,u[k][j][i][15]*user->p0,u[k][j][i][16]*user->B0,u[k][j][i][17]*user->B0,u[k][j][i][18]*user->B0);
          PetscPrintf(PETSC_COMM_SELF,"@[%3d][%3d][%3d] ve=[%+9.3e %+9.3e %+9.3e]; J=[%+9.3e %+9.3e %+9.3e]; E=[%+9.3e %+9.3e %+9.3e]\n",i,j,k,v[k][j][i][0]*user->v0,v[k][j][i][1]*user->v0,v[k][j][i][2]*user->v0,v[k][j][i][3]*qe*user->n0*user->v0,v[k][j][i][4]*qe*user->n0*user->v0,v[k][j][i][5]*qe*user->n0*user->v0,v[k][j][i][6]*user->v0*user->B0,v[k][j][i][7]*user->v0*user->B0,v[k][j][i][8]*user->v0*user->B0);
          exit(1);
        }
        */
      }
    }
  }

  m = MPI_Allreduce(MPI_IN_PLACE,&ve_max[0]   ,3 ,MPIU_REAL,MPIU_MAX,PETSC_COMM_WORLD);
  m = MPI_Allreduce(MPI_IN_PLACE,&vi_max[0][0],9 ,MPIU_REAL,MPIU_MAX,PETSC_COMM_WORLD);
  m = MPI_Allreduce(MPI_IN_PLACE,&vf_max.g[0] ,51,MPIU_REAL,MPIU_MAX,PETSC_COMM_WORLD);
  m = MPI_Allreduce(MPI_IN_PLACE,&dtCFL       ,1 ,MPIU_REAL,MPIU_MIN,PETSC_COMM_WORLD);
  //ierr = PetscPrintf(PETSC_COMM_SELF,"@%d ta = %+14.7e ts_n = %+14.7e ts_c = %+14.7e \n",rank, tA*tau,ts_n*tau,ts_c*tau);CHKERRQ(ierr);
  m = MPI_Allreduce(MPI_IN_PLACE,&tA          ,1 ,MPIU_REAL,MPIU_SUM,PETSC_COMM_WORLD);
  m = MPI_Allreduce(MPI_IN_PLACE,&ts_n        ,1 ,MPIU_REAL,MPIU_SUM,PETSC_COMM_WORLD);
  m = MPI_Allreduce(MPI_IN_PLACE,&ts_c        ,1 ,MPIU_REAL,MPIU_SUM,PETSC_COMM_WORLD);
  //ierr = PetscPrintf(PETSC_COMM_WORLD,"ta = %+14.7e ts_n = %+14.7e ts_c = %+14.7e \n",tA*tau,ts_n*tau,ts_c*tau);CHKERRQ(ierr);
  //exit(12);

  //ierr = PetscPrintf(PETSC_COMM_WORLD,"\t==> t=%12.5e: dt.CFL=%12.5e:\nv_max[O2p]    =[%+12.3e %+12.3e %+12.3e] v_max[CO2p]    =[%+12.3e %+12.3e %+12.3e] v_max[Op]    =[%+12.3e %+12.3e %+12.3e] v_max[e]    =[%+12.3e %+12.3e %+12.3e]\nff_max.g      =[%+12.3e %+12.3e %+12.3e]\nff_max.dp[O2p]=[%+12.3e %+12.3e %+12.3e] ff_max.dp[CO2p]=[%+12.3e %+12.3e %+12.3e] ff_max.dp[Op]=[%+12.3e %+12.3e %+12.3e] ff_max.dp[e]=[%+12.3e %+12.3e %+12.3e]\nff_max.E[O2p] =[%+12.3e %+12.3e %+12.3e] ff_max.E[CO2p] =[%+12.3e %+12.3e %+12.3e] ff_max.E[Op] =[%+12.3e %+12.3e %+12.3e] ff_max.E[e] =[%+12.3e %+12.3e %+12.3e]\nff_max.B[O2p] =[%+12.3e %+12.3e %+12.3e] ff_max.B[CO2p] =[%+12.3e %+12.3e %+12.3e] ff_max.B[Op] =[%+12.3e %+12.3e %+12.3e] ff_max.B[e] =[%+12.3e %+12.3e %+12.3e]\nff_max.c[O2p] =[%+12.3e %+12.3e %+12.3e] ff_max.c[CO2p] =[%+12.3e %+12.3e %+12.3e] ff_max.c[Op] =[%+12.3e %+12.3e %+12.3e] ff_max.c[e] =[%+12.3e %+12.3e %+12.3e]\n",ftime*tau,dt*user->tau,vi_max[O2p][0]*v0,vi_max[O2p][1]*v0,vi_max[O2p][2]*v0,vi_max[CO2p][0]*v0,vi_max[CO2p][1]*v0,vi_max[CO2p][2]*v0,vi_max[Op][0]*v0,vi_max[Op][1]*v0,vi_max[Op][2]*v0,ve_max[0]*v0,ve_max[1]*v0,ve_max[2]*v0,ff_max.g[0]*v0*user->B0,ff_max.g[1]*v0*user->B0,ff_max.g[2]*v0*user->B0,ff_max.dp[O2p][0]*v0*user->B0,ff_max.dp[O2p][1]*v0*user->B0,ff_max.dp[O2p][2]*v0*user->B0,ff_max.dp[CO2p][0]*v0*user->B0,ff_max.dp[CO2p][1]*v0*user->B0,ff_max.dp[CO2p][2]*v0*user->B0,ff_max.dp[Op][0]*v0*user->B0,ff_max.dp[Op][1]*v0*user->B0,ff_max.dp[Op][2]*v0*user->B0,ff_max.dp[e][0]*v0*user->B0,ff_max.dp[e][1]*v0*user->B0,ff_max.dp[e][2]*v0*user->B0,ff_max.E[O2p][0]*v0*user->B0,ff_max.E[O2p][1]*v0*user->B0,ff_max.E[O2p][2]*v0*user->B0,ff_max.E[CO2p][0]*v0*user->B0,ff_max.E[CO2p][1]*v0*user->B0,ff_max.E[CO2p][2]*v0*user->B0,ff_max.E[Op][0]*v0*user->B0,ff_max.E[Op][1]*v0*user->B0,ff_max.E[Op][2]*v0*user->B0,ff_max.E[e][0]*v0*user->B0,ff_max.E[e][1]*v0*user->B0,ff_max.E[e][2]*v0*user->B0,ff_max.B[O2p][0]*v0*user->B0,ff_max.B[O2p][1]*v0*user->B0,ff_max.B[O2p][2]*v0*user->B0,ff_max.B[CO2p][0]*v0*user->B0,ff_max.B[CO2p][1]*v0*user->B0,ff_max.B[CO2p][2]*v0*user->B0,ff_max.B[Op][0]*v0*user->B0,ff_max.B[Op][1]*v0*user->B0,ff_max.B[Op][2]*v0*user->B0,ff_max.B[e][0]*v0*user->B0,ff_max.B[e][1]*v0*user->B0,ff_max.B[e][2]*v0*user->B0,ff_max.c[O2p][0]*v0*user->B0,ff_max.c[O2p][1]*v0*user->B0,ff_max.c[O2p][2]*v0*user->B0,ff_max.c[CO2p][0]*v0*user->B0,ff_max.c[CO2p][1]*v0*user->B0,ff_max.c[CO2p][2]*v0*user->B0,ff_max.c[Op][0]*v0*user->B0,ff_max.c[Op][1]*v0*user->B0,ff_max.c[Op][2]*v0*user->B0,ff_max.c[e][0]*v0*user->B0,ff_max.c[e][1]*v0*user->B0,ff_max.c[e][2]*v0*user->B0); CHKERRQ(ierr);
  
  /*
   * BELOW IS THE LIST OF DIAGNOSTIC TO ACTIVATE
   */
  /*
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\t==> t=%12.5e: dt.CFL=%12.5e MAX VALUES: \ntA = %+14.7e ts_n = %+14.7e ts_c = %+14.7e\nv.O2p    =[%+14.7e %+14.7e %+14.7e]; v.CO2p    =[%+14.7e %+14.7e %+14.7e]; v.Op    =[%+14.7e %+14.7e %+14.7e]; v.e    =[%+14.7e %+14.7e %+14.7e]\n", 
              ftime*tau,dt*tau,
              tA                 *tau   ,ts_n               *tau   ,ts_c               *tau   ,
              vi_max[O2p][0]     *v0    ,vi_max[O2p][1]     *v0    ,vi_max[O2p][2]     *v0    ,
              vi_max[CO2p][0]    *v0    ,vi_max[CO2p][1]    *v0    ,vi_max[CO2p][2]    *v0    ,
              vi_max[Op][0]      *v0    ,vi_max[Op][1]      *v0    ,vi_max[Op][2]      *v0    ,
              ve_max[0]          *v0    ,ve_max[1]          *v0    ,ve_max[2]          *v0    );CHKERRQ(ierr);
   */
  /*
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\t==> t=%12.5e: dt.CFL=%12.5e MAX VALUES: \nv.O2p    =[%+14.7e %+14.7e %+14.7e]; v.CO2p    =[%+14.7e %+14.7e %+14.7e]; v.Op    =[%+14.7e %+14.7e %+14.7e]; v.e    =[%+14.7e %+14.7e %+14.7e]\na.g      =[%+14.7e %+14.7e %+14.7e]\na.dp.O2p =[%+14.7e %+14.7e %+14.7e]; a.dp.CO2p =[%+14.7e %+14.7e %+14.7e]; a.dp.Op =[%+14.7e %+14.7e %+14.7e]; a.dp.e =[%+14.7e %+14.7e %+14.7e]\na.E.O2p  =[%+14.7e %+14.7e %+14.7e]; a.E.CO2p  =[%+14.7e %+14.7e %+14.7e]; a.E.Op  =[%+14.7e %+14.7e %+14.7e]; a.E.e  =[%+14.7e %+14.7e %+14.7e]\na.B.O2p  =[%+14.7e %+14.7e %+14.7e]; a.B.CO2p  =[%+14.7e %+14.7e %+14.7e]; a.B.Op  =[%+14.7e %+14.7e %+14.7e]; a.B.e  =[%+14.7e %+14.7e %+14.7e]\na.c.O2p  =[%+14.7e %+14.7e %+14.7e]; a.c.CO2p  =[%+14.7e %+14.7e %+14.7e]; a.c.Op  =[%+14.7e %+14.7e %+14.7e]; a.c.e  =[%+14.7e %+14.7e %+14.7e]\na.adv.O2p=[%+14.7e %+14.7e %+14.7e]; a.adv.CO2p=[%+14.7e %+14.7e %+14.7e]; a.adv.Op=[%+14.7e %+14.7e %+14.7e]; a.adv.e=[%+14.7e %+14.7e %+14.7e]\n", 
              ftime              *tau   ,dt                 *tau   ,
              vi_max[O2p][0]     *v0    ,vi_max[O2p][1]     *v0    ,vi_max[O2p][2]     *v0    ,
              vi_max[CO2p][0]    *v0    ,vi_max[CO2p][1]    *v0    ,vi_max[CO2p][2]    *v0    ,
              vi_max[Op][0]      *v0    ,vi_max[Op][1]      *v0    ,vi_max[Op][2]      *v0    ,
              ve_max[0]          *v0    ,ve_max[1]          *v0    ,ve_max[2]          *v0    ,
              ff_max.g[0]        *v0/tau,ff_max.g[1]        *v0/tau,ff_max.g[2]        *v0/tau,
              ff_max.dp[O2p][0]  *v0/tau,ff_max.dp[O2p][1]  *v0/tau,ff_max.dp[O2p][2]  *v0/tau,
              ff_max.dp[CO2p][0] *v0/tau,ff_max.dp[CO2p][1] *v0/tau,ff_max.dp[CO2p][2] *v0/tau,
              ff_max.dp[Op][0]   *v0/tau,ff_max.dp[Op][1]   *v0/tau,ff_max.dp[Op][2]   *v0/tau,
              ff_max.dp[e][0]    *v0/tau,ff_max.dp[e][1]    *v0/tau,ff_max.dp[e][2]    *v0/tau,
              ff_max.E[O2p][0]   *v0/tau,ff_max.E[O2p][1]   *v0/tau,ff_max.E[O2p][2]   *v0/tau,
              ff_max.E[CO2p][0]  *v0/tau,ff_max.E[CO2p][1]  *v0/tau,ff_max.E[CO2p][2]  *v0/tau,
              ff_max.E[Op][0]    *v0/tau,ff_max.E[Op][1]    *v0/tau,ff_max.E[Op][2]    *v0/tau,
              ff_max.E[e][0]     *v0/tau,ff_max.E[e][1]     *v0/tau,ff_max.E[e][2]     *v0/tau,
              ff_max.B[O2p][0]   *v0/tau,ff_max.B[O2p][1]   *v0/tau,ff_max.B[O2p][2]   *v0/tau,
              ff_max.B[CO2p][0]  *v0/tau,ff_max.B[CO2p][1]  *v0/tau,ff_max.B[CO2p][2]  *v0/tau,
              ff_max.B[Op][0]    *v0/tau,ff_max.B[Op][1]    *v0/tau,ff_max.B[Op][2]    *v0/tau,
              ff_max.c[O2p][0]   *v0/tau,ff_max.c[O2p][1]   *v0/tau,ff_max.c[O2p][2]   *v0/tau,
              ff_max.c[CO2p][0]  *v0/tau,ff_max.c[CO2p][1]  *v0/tau,ff_max.c[CO2p][2]  *v0/tau,
              ff_max.c[Op][0]    *v0/tau,ff_max.c[Op][1]    *v0/tau,ff_max.c[Op][2]    *v0/tau,
              ff_max.c[e][0]     *v0/tau,ff_max.c[e][1]     *v0/tau,ff_max.c[e][2]     *v0/tau,
              ff_max.adv[O2p][0] *v0/tau,ff_max.adv[O2p][1] *v0/tau,ff_max.adv[O2p][2] *v0/tau,
              ff_max.adv[CO2p][0]*v0/tau,ff_max.adv[CO2p][1]*v0/tau,ff_max.adv[CO2p][2]*v0/tau,
              ff_max.adv[Op][0]  *v0/tau,ff_max.adv[Op][1]  *v0/tau,ff_max.adv[Op][2]  *v0/tau,
              ff_max.adv[e][0]   *v0/tau,ff_max.adv[e][1]   *v0/tau,ff_max.adv[e][2]   *v0/tau
                ); CHKERRQ(ierr);
    */
  /*
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\t==> MAX DRIFT VALUES: \nv.O2p    =[%+14.7e %+14.7e %+14.7e]; v.CO2p    =[%+14.7e %+14.7e %+14.7e]; v.Op    =[%+14.7e %+14.7e %+14.7e]; v.e    =[%+14.7e %+14.7e %+14.7e]\nv.g      =[%+14.7e %+14.7e %+14.7e]\nv.dp.O2p =[%+14.7e %+14.7e %+14.7e]; v.dp.CO2p =[%+14.7e %+14.7e %+14.7e]; v.dp.Op =[%+14.7e %+14.7e %+14.7e]; v.dp.e =[%+14.7e %+14.7e %+14.7e]\nv.E.O2p  =[%+14.7e %+14.7e %+14.7e]; v.E.CO2p  =[%+14.7e %+14.7e %+14.7e]; v.E.Op  =[%+14.7e %+14.7e %+14.7e]; v.E.e  =[%+14.7e %+14.7e %+14.7e]\nv.B.O2p  =[%+14.7e %+14.7e %+14.7e]; v.B.CO2p  =[%+14.7e %+14.7e %+14.7e]; v.B.Op  =[%+14.7e %+14.7e %+14.7e]; v.B.e  =[%+14.7e %+14.7e %+14.7e]\nv.c.O2p  =[%+14.7e %+14.7e %+14.7e]; v.c.CO2p  =[%+14.7e %+14.7e %+14.7e]; v.c.Op  =[%+14.7e %+14.7e %+14.7e]; v.c.e  =[%+14.7e %+14.7e %+14.7e]\nv.adv.O2p=[%+14.7e %+14.7e %+14.7e]; v.adv.CO2p=[%+14.7e %+14.7e %+14.7e]; v.adv.Op=[%+14.7e %+14.7e %+14.7e]; v.adv.e=[%+14.7e %+14.7e %+14.7e]\n", 
              vi_max[O2p][0]     *v0,vi_max[O2p][1]     *v0,vi_max[O2p][2]     *v0,
              vi_max[CO2p][0]    *v0,vi_max[CO2p][1]    *v0,vi_max[CO2p][2]    *v0,
              vi_max[Op][0]      *v0,vi_max[Op][1]      *v0,vi_max[Op][2]      *v0,
              ve_max[0]          *v0,ve_max[1]          *v0,ve_max[2]          *v0,
              vf_max.g[0]        *v0,vf_max.g[1]        *v0,vf_max.g[2]        *v0,
              vf_max.dp[O2p][0]  *v0,vf_max.dp[O2p][1]  *v0,vf_max.dp[O2p][2]  *v0,
              vf_max.dp[CO2p][0] *v0,vf_max.dp[CO2p][1] *v0,vf_max.dp[CO2p][2] *v0,
              vf_max.dp[Op][0]   *v0,vf_max.dp[Op][1]   *v0,vf_max.dp[Op][2]   *v0,
              vf_max.dp[e][0]    *v0,vf_max.dp[e][1]    *v0,vf_max.dp[e][2]    *v0,
              vf_max.E[O2p][0]   *v0,vf_max.E[O2p][1]   *v0,vf_max.E[O2p][2]   *v0,
              vf_max.E[CO2p][0]  *v0,vf_max.E[CO2p][1]  *v0,vf_max.E[CO2p][2]  *v0,
              vf_max.E[Op][0]    *v0,vf_max.E[Op][1]    *v0,vf_max.E[Op][2]    *v0,
              vf_max.E[e][0]     *v0,vf_max.E[e][1]     *v0,vf_max.E[e][2]     *v0,
              vf_max.B[O2p][0]   *v0,vf_max.B[O2p][1]   *v0,vf_max.B[O2p][2]   *v0,
              vf_max.B[CO2p][0]  *v0,vf_max.B[CO2p][1]  *v0,vf_max.B[CO2p][2]  *v0,
              vf_max.B[Op][0]    *v0,vf_max.B[Op][1]    *v0,vf_max.B[Op][2]    *v0,
              vf_max.c[O2p][0]   *v0,vf_max.c[O2p][1]   *v0,vf_max.c[O2p][2]   *v0,
              vf_max.c[CO2p][0]  *v0,vf_max.c[CO2p][1]  *v0,vf_max.c[CO2p][2]  *v0,
              vf_max.c[Op][0]    *v0,vf_max.c[Op][1]    *v0,vf_max.c[Op][2]    *v0,
              vf_max.c[e][0]     *v0,vf_max.c[e][1]     *v0,vf_max.c[e][2]     *v0,
              vf_max.adv[O2p][0] *v0,vf_max.adv[O2p][1] *v0,vf_max.adv[O2p][2] *v0,
              vf_max.adv[CO2p][0]*v0,vf_max.adv[CO2p][1]*v0,vf_max.adv[CO2p][2]*v0,
              vf_max.adv[Op][0]  *v0,vf_max.adv[Op][1]  *v0,vf_max.adv[Op][2]  *v0,
              vf_max.adv[e][0]   *v0,vf_max.adv[e][1]   *v0,vf_max.adv[e][2]   *v0
                ); CHKERRQ(ierr);
*/
  // Adapt timestep to satisfy CFL //
//user->dt = 0.01*dtCFL;

  /*
  if(rank==0)
    PetscPrintf(PETSC_COMM_WORLD,"step = %d\n",user->cnt);
  if (rank == 7 && user->cnt==20) { 
    k=mz-1;
    j=my-1;
    i=mx-1; 
    PetscPrintf(PETSC_COMM_SELF,"c.nO2p=%+14.7e; c.nCO2p=%+14.7e; c.nOp=%+14.7e; c.vO2p=[%+14.7e %+14.7e %+14.7e]; c.vCO2p=[%+14.7e %+14.7e %+14.7e]; c.vOp=[%+14.7e %+14.7e %+14.7e]; c.pO2p=%+14.7e; c.pCO2p=%+14.7e; c.pOp=%+14.7e; c.pe=%+14.7e; c.B=[%+14.7e %+14.7e %+14.7e];\n",i,j,k,u[k][j][i][0]*user->n0,u[k][j][i][1]*user->n0,u[k][j][i][2]*user->n0,u[k][j][i][3]*user->v0,u[k][j][i][4]*user->v0,u[k][j][i][5]*user->v0,u[k][j][i][6]*user->v0,u[k][j][i][7]*user->v0,u[k][j][i][8]*user->v0,u[k][j][i][9]*user->v0,u[k][j][i][10]*user->v0,u[k][j][i][11]*user->v0,u[k][j][i][12]*user->p0,u[k][j][i][13]*user->p0,u[k][j][i][14]*user->p0,u[k][j][i][15]*user->p0,u[k][j][i][16]*user->B0,u[k][j][i][17]*user->B0,u[k][j][i][18]*user->B0);
    PetscPrintf(PETSC_COMM_SELF,"c.ve=[%+14.7e %+14.7e %+14.7e]; c.J=[%+14.7e %+14.7e %+14.7e]; c.E=[%+14.7e %+14.7e %+14.7e];\n",i,j,k,v[k][j][i][0]*user->v0,v[k][j][i][1]*user->v0,v[k][j][i][2]*user->v0,v[k][j][i][3]*qe*user->n0*user->v0,v[k][j][i][4]*qe*user->n0*user->v0,v[k][j][i][5]*qe*user->n0*user->v0,v[k][j][i][6]*user->v0*user->B0,v[k][j][i][7]*user->v0*user->B0,v[k][j][i][8]*user->v0*user->B0);
    PetscPrintf(PETSC_COMM_SELF,"u.nO2p=%+14.7e; u.nCO2p=%+14.7e; u.nOp=%+14.7e; u.vO2p=[%+14.7e %+14.7e %+14.7e]; u.vCO2p=[%+14.7e %+14.7e %+14.7e]; u.vOp=[%+14.7e %+14.7e %+14.7e]; u.pO2p=%+14.7e; u.pCO2p=%+14.7e; u.pOp=%+14.7e; u.pe=%+14.7e; u.B=[%+14.7e %+14.7e %+14.7e];\n",i,j,k,u[k+1][j][i][0]*user->n0,u[k+1][j][i][1]*user->n0,u[k+1][j][i][2]*user->n0,u[k+1][j][i][3]*user->v0,u[k+1][j][i][4]*user->v0,u[k+1][j][i][5]*user->v0,u[k+1][j][i][6]*user->v0,u[k+1][j][i][7]*user->v0,u[k+1][j][i][8]*user->v0,u[k+1][j][i][9]*user->v0,u[k+1][j][i][10]*user->v0,u[k+1][j][i][11]*user->v0,u[k+1][j][i][12]*user->p0,u[k+1][j][i][13]*user->p0,u[k+1][j][i][14]*user->p0,u[k+1][j][i][15]*user->p0,u[k+1][j][i][16]*user->B0,u[k+1][j][i][17]*user->B0,u[k+1][j][i][18]*user->B0);
    PetscPrintf(PETSC_COMM_SELF,"u.ve=[%+14.7e %+14.7e %+14.7e]; u.J=[%+14.7e %+14.7e %+14.7e]; u.E=[%+14.7e %+14.7e %+14.7e];\n",i,j,k,v[k+1][j][i][0]*user->v0,v[k+1][j][i][1]*user->v0,v[k+1][j][i][2]*user->v0,v[k+1][j][i][3]*qe*user->n0*user->v0,v[k+1][j][i][4]*qe*user->n0*user->v0,v[k+1][j][i][5]*qe*user->n0*user->v0,v[k+1][j][i][6]*user->v0*user->B0,v[k+1][j][i][7]*user->v0*user->B0,v[k+1][j][i][8]*user->v0*user->B0);
    PetscPrintf(PETSC_COMM_SELF,"d.nO2p=%+14.7e; d.nCO2p=%+14.7e; d.nOp=%+14.7e; d.vO2p=[%+14.7e %+14.7e %+14.7e]; d.vCO2p=[%+14.7e %+14.7e %+14.7e]; d.vOp=[%+14.7e %+14.7e %+14.7e]; d.pO2p=%+14.7e; d.pCO2p=%+14.7e; d.pOp=%+14.7e; d.pe=%+14.7e; d.B=[%+14.7e %+14.7e %+14.7e];\n",i,j,k,u[k-1][j][i][0]*user->n0,u[k-1][j][i][1]*user->n0,u[k-1][j][i][2]*user->n0,u[k-1][j][i][3]*user->v0,u[k-1][j][i][4]*user->v0,u[k-1][j][i][5]*user->v0,u[k-1][j][i][6]*user->v0,u[k-1][j][i][7]*user->v0,u[k-1][j][i][8]*user->v0,u[k-1][j][i][9]*user->v0,u[k-1][j][i][10]*user->v0,u[k-1][j][i][11]*user->v0,u[k-1][j][i][12]*user->p0,u[k-1][j][i][13]*user->p0,u[k-1][j][i][14]*user->p0,u[k-1][j][i][15]*user->p0,u[k-1][j][i][16]*user->B0,u[k-1][j][i][17]*user->B0,u[k-1][j][i][18]*user->B0);
    PetscPrintf(PETSC_COMM_SELF,"d.ve=[%+14.7e %+14.7e %+14.7e]; d.J=[%+14.7e %+14.7e %+14.7e]; d.E=[%+14.7e %+14.7e %+14.7e];\n",i,j,k,v[k-1][j][i][0]*user->v0,v[k-1][j][i][1]*user->v0,v[k-1][j][i][2]*user->v0,v[k-1][j][i][3]*qe*user->n0*user->v0,v[k-1][j][i][4]*qe*user->n0*user->v0,v[k-1][j][i][5]*qe*user->n0*user->v0,v[k-1][j][i][6]*user->v0*user->B0,v[k-1][j][i][7]*user->v0*user->B0,v[k-1][j][i][8]*user->v0*user->B0);
    PetscPrintf(PETSC_COMM_SELF,"f.nO2p=%+14.7e; f.nCO2p=%+14.7e; f.nOp=%+14.7e; f.vO2p=[%+14.7e %+14.7e %+14.7e]; f.vCO2p=[%+14.7e %+14.7e %+14.7e]; f.vOp=[%+14.7e %+14.7e %+14.7e]; f.pO2p=%+14.7e; f.pCO2p=%+14.7e; f.pOp=%+14.7e; f.pe=%+14.7e; f.B=[%+14.7e %+14.7e %+14.7e];\n",i,j,k,u[k][j][i+1][0]*user->n0,u[k][j][i+1][1]*user->n0,u[k][j][i+1][2]*user->n0,u[k][j][i+1][3]*user->v0,u[k][j][i+1][4]*user->v0,u[k][j][i+1][5]*user->v0,u[k][j][i+1][6]*user->v0,u[k][j][i+1][7]*user->v0,u[k][j][i+1][8]*user->v0,u[k][j][i+1][9]*user->v0,u[k][j][i+1][10]*user->v0,u[k][j][i+1][11]*user->v0,u[k][j][i+1][12]*user->p0,u[k][j][i+1][13]*user->p0,u[k][j][i+1][14]*user->p0,u[k][j][i+1][15]*user->p0,u[k][j][i+1][16]*user->B0,u[k][j][i+1][17]*user->B0,u[k][j][i+1][18]*user->B0);
    PetscPrintf(PETSC_COMM_SELF,"f.ve=[%+14.7e %+14.7e %+14.7e]; f.J=[%+14.7e %+14.7e %+14.7e]; f.E=[%+14.7e %+14.7e %+14.7e];\n",i,j,k,v[k][j][i+1][0]*user->v0,v[k][j][i+1][1]*user->v0,v[k][j][i+1][2]*user->v0,v[k][j][i+1][3]*qe*user->n0*user->v0,v[k][j][i+1][4]*qe*user->n0*user->v0,v[k][j][i+1][5]*qe*user->n0*user->v0,v[k][j][i+1][6]*user->v0*user->B0,v[k][j][i+1][7]*user->v0*user->B0,v[k][j][i+1][8]*user->v0*user->B0);
    PetscPrintf(PETSC_COMM_SELF,"b.nO2p=%+14.7e; b.nCO2p=%+14.7e; b.nOp=%+14.7e; b.vO2p=[%+14.7e %+14.7e %+14.7e]; b.vCO2p=[%+14.7e %+14.7e %+14.7e]; b.vOp=[%+14.7e %+14.7e %+14.7e]; b.pO2p=%+14.7e; b.pCO2p=%+14.7e; b.pOp=%+14.7e; b.pe=%+14.7e; b.B=[%+14.7e %+14.7e %+14.7e];\n",i,j,k,u[k][j][i-1][0]*user->n0,u[k][j][i-1][1]*user->n0,u[k][j][i-1][2]*user->n0,u[k][j][i-1][3]*user->v0,u[k][j][i-1][4]*user->v0,u[k][j][i-1][5]*user->v0,u[k][j][i-1][6]*user->v0,u[k][j][i-1][7]*user->v0,u[k][j][i-1][8]*user->v0,u[k][j][i-1][9]*user->v0,u[k][j][i-1][10]*user->v0,u[k][j][i-1][11]*user->v0,u[k][j][i-1][12]*user->p0,u[k][j][i-1][13]*user->p0,u[k][j][i-1][14]*user->p0,u[k][j][i-1][15]*user->p0,u[k][j][i-1][16]*user->B0,u[k][j][i-1][17]*user->B0,u[k][j][i-1][18]*user->B0);
    PetscPrintf(PETSC_COMM_SELF,"b.ve=[%+14.7e %+14.7e %+14.7e]; b.J=[%+14.7e %+14.7e %+14.7e]; b.E=[%+14.7e %+14.7e %+14.7e];\n",i,j,k,v[k][j][i-1][0]*user->v0,v[k][j][i-1][1]*user->v0,v[k][j][i-1][2]*user->v0,v[k][j][i-1][3]*qe*user->n0*user->v0,v[k][j][i-1][4]*qe*user->n0*user->v0,v[k][j][i-1][5]*qe*user->n0*user->v0,v[k][j][i-1][6]*user->v0*user->B0,v[k][j][i-1][7]*user->v0*user->B0,v[k][j][i-1][8]*user->v0*user->B0);
    PetscPrintf(PETSC_COMM_SELF,"r.nO2p=%+14.7e; r.nCO2p=%+14.7e; r.nOp=%+14.7e; r.vO2p=[%+14.7e %+14.7e %+14.7e]; r.vCO2p=[%+14.7e %+14.7e %+14.7e]; r.vOp=[%+14.7e %+14.7e %+14.7e]; r.pO2p=%+14.7e; r.pCO2p=%+14.7e; r.pOp=%+14.7e; r.pe=%+14.7e; r.B=[%+14.7e %+14.7e %+14.7e];\n",i,j,k,u[k][j+1][i][0]*user->n0,u[k][j+1][i][1]*user->n0,u[k][j+1][i][2]*user->n0,u[k][j+1][i][3]*user->v0,u[k][j+1][i][4]*user->v0,u[k][j+1][i][5]*user->v0,u[k][j+1][i][6]*user->v0,u[k][j+1][i][7]*user->v0,u[k][j+1][i][8]*user->v0,u[k][j+1][i][9]*user->v0,u[k][j+1][i][10]*user->v0,u[k][j+1][i][11]*user->v0,u[k][j+1][i][12]*user->p0,u[k][j+1][i][13]*user->p0,u[k][j+1][i][14]*user->p0,u[k][j+1][i][15]*user->p0,u[k][j+1][i][16]*user->B0,u[k][j+1][i][17]*user->B0,u[k][j+1][i][18]*user->B0);
    PetscPrintf(PETSC_COMM_SELF,"r.ve=[%+14.7e %+14.7e %+14.7e]; r.J=[%+14.7e %+14.7e %+14.7e]; r.E=[%+14.7e %+14.7e %+14.7e];\n",i,j,k,v[k][j+1][i][0]*user->v0,v[k][j+1][i][1]*user->v0,v[k][j+1][i][2]*user->v0,v[k][j+1][i][3]*qe*user->n0*user->v0,v[k][j+1][i][4]*qe*user->n0*user->v0,v[k][j+1][i][5]*qe*user->n0*user->v0,v[k][j+1][i][6]*user->v0*user->B0,v[k][j+1][i][7]*user->v0*user->B0,v[k][j+1][i][8]*user->v0*user->B0);
    PetscPrintf(PETSC_COMM_SELF,"l.nO2p=%+14.7e; l.nCO2p=%+14.7e; l.nOp=%+14.7e; l.vO2p=[%+14.7e %+14.7e %+14.7e]; l.vCO2p=[%+14.7e %+14.7e %+14.7e]; l.vOp=[%+14.7e %+14.7e %+14.7e]; l.pO2p=%+14.7e; l.pCO2p=%+14.7e; l.pOp=%+14.7e; l.pe=%+14.7e; l.B=[%+14.7e %+14.7e %+14.7e];\n",i,j,k,u[k][j-1][i][0]*user->n0,u[k][j-1][i][1]*user->n0,u[k][j-1][i][2]*user->n0,u[k][j-1][i][3]*user->v0,u[k][j-1][i][4]*user->v0,u[k][j-1][i][5]*user->v0,u[k][j-1][i][6]*user->v0,u[k][j-1][i][7]*user->v0,u[k][j-1][i][8]*user->v0,u[k][j-1][i][9]*user->v0,u[k][j-1][i][10]*user->v0,u[k][j-1][i][11]*user->v0,u[k][j-1][i][12]*user->p0,u[k][j-1][i][13]*user->p0,u[k][j-1][i][14]*user->p0,u[k][j-1][i][15]*user->p0,u[k][j-1][i][16]*user->B0,u[k][j-1][i][17]*user->B0,u[k][j-1][i][18]*user->B0);
    PetscPrintf(PETSC_COMM_SELF,"l.ve=[%+14.7e %+14.7e %+14.7e]; l.J=[%+14.7e %+14.7e %+14.7e]; l.E=[%+14.7e %+14.7e %+14.7e];\n",i,j,k,v[k][j-1][i][0]*user->v0,v[k][j-1][i][1]*user->v0,v[k][j-1][i][2]*user->v0,v[k][j-1][i][3]*qe*user->n0*user->v0,v[k][j-1][i][4]*qe*user->n0*user->v0,v[k][j-1][i][5]*qe*user->n0*user->v0,v[k][j-1][i][6]*user->v0*user->B0,v[k][j-1][i][7]*user->v0*user->B0,v[k][j-1][i][8]*user->v0*user->B0);
    PetscPrintf(PETSC_COMM_SELF,"c.F.nO2p=%+14.7e; c.F.nCO2p=%+14.7e; c.F.nOp=%+14.7e; c.F.vO2p=[%+14.7e %+14.7e %+14.7e]; c.F.vCO2p=[%+14.7e %+14.7e %+14.7e]; c.F.vOp=[%+14.7e %+14.7e %+14.7e]; c.F.pO2p=%+14.7e; c.F.pCO2p=%+14.7e; c.F.pOp=%+14.7e; c.F.pe=%+14.7e; c.F.B=[%+14.7e %+14.7e %+14.7e];\n",i,j,k,f[k][j][i][0]*user->n0/user->tau,f[k][j][i][1]*user->n0/user->tau,f[k][j][i][2]*user->n0/user->tau,f[k][j][i][3]*user->v0/user->tau,f[k][j][i][4]*user->v0/user->tau,f[k][j][i][5]*user->v0/user->tau,f[k][j][i][6]*user->v0/user->tau,f[k][j][i][7]*user->v0/user->tau,f[k][j][i][8]*user->v0/user->tau,f[k][j][i][9]*user->v0/user->tau,f[k][j][i][10]*user->v0/user->tau,f[k][j][i][11]*user->v0/user->tau,f[k][j][i][12]*user->p0/user->tau,f[k][j][i][13]*user->p0/user->tau,f[k][j][i][14]*user->p0/user->tau,f[k][j][i][15]*user->p0/user->tau,f[k][j][i][16]*user->B0/user->tau,f[k][j][i][17]*user->B0/user->tau,f[k][j][i][18]*user->B0);
    exit(12);
  }
  */

  // Restore vectors //
  ierr = DMRestoreGlobalVector(db,&V);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(db,localV,&v);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(da,localU,&u);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(da,F,&f);CHKERRQ(ierr);
  
  ierr = DMRestoreLocalVector(db,&localV);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(da,localU,INSERT_VALUES,U);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(da,localU,INSERT_VALUES,U);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(da,&localU);CHKERRQ(ierr);
  PetscFunctionReturn(0); 
} 

