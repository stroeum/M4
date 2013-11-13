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
  PetscReal      L=user->L,Lz;
  PetscReal      Xmin=user->outXmin, Ymin=user->outYmin, Zmin=user->outZmin;
  PetscReal      inZmax=user->inZmax;
//PetscReal      /*me = user->me,*/ mi[3] = {user->mi[O2p], user->mi[CO2p], user->mi[Op]};
  PetscReal      n0=user->n0, v0=user->v0, p0=user->p0, B0=user->B0;
  PetscReal      Bo=user->B[0], a=user->B[1], b=user->B[2], c=user->B[3];
//PetscReal      un[3]={user->un[0],user->un[1],user->un[2]}; // neutral wind
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
  PetscReal      nio[3],neo/*,vio[3]*/,pio[3],peo;
  PetscInt       rank;
  PetscViewer    fViewer;
  char           fName[PETSC_MAX_PATH_LEN];
  PetscBool      flag;

  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"Form Initial Conditions: ...\n");
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  Lz = user->outZmax-user->inZmin;
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
    // Map U //
    ierr = DMDAVecGetArrayDOF(da,U,&u);CHKERRQ(ierr);

    // Get local grid boundaries //
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

          /*
          Te = 2000;
          Ti = 2000;
          neo = 1e11;
          nio[O2p]  = .71*neo;
          nio[CO2p] = .25*neo;
          nio[Op]   = .04*neo;
          //if(!(neo>0)) neo = eps.n*n0;
          */

          for (l=0; l<3; l++) {
            if (!(nio[l] > 0)) nio[l] = eps.n*n0;
            //vio[l] = 100.0;
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

              //if (k==13) if (j==24) PetscPrintf(PETSC_COMM_SELF,"@%d [%+12.6f %+12.6f %+12.6f] %f\n",rank,X,Y,Z,ui[m]);
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
            if (Btype==2 && ui[2]==0) {
              //PetscPrintf(PETSC_COMM_WORLD,"[%i %i %i], It worked... so far.\n",i,j,k);
              //SETERRQ(PETSC_COMM_WORLD,62,"For a horizontal field, an initial vertical velocity is preferred.");
            }
            //vmax  = PetscMax(vmax,vio[l]/v0);
            u[k][j][i][d.pi[l]] = pio[l]/p0;
          }
          //vmax = PetscMax(vmax,veo/v0);
          u[k][j][i][d.pe] = peo/p0;
        }
      }
    }

    ierr = DMDAVecRestoreArrayDOF(da,U,&u);CHKERRQ(ierr);
  }

  //PetscPrintf(PETSC_COMM_WORLD,"ARE WE THERE YET?\n");
  // Compute initial BC for u //
  //ierr = FormBCu(U,user);CHKERRQ(ierr);

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


  for (k=zs; k<zs+zm; k++) {
    Z = Zmin + z[k]*L;
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
        if (bcType==11) { // Combination of DiRichlet and Neumann
          for (l=0; l< 3; l++) u[k][j][id.x[0]][l] = Interpolate(user->RefPart, l ,Z, lin_flat)*Interpolate(user->RefProf, 7 ,Z, lin_exp)/n0; 
          for (l=3; l<19; l++) u[k][j][id.x[0]][l] = u[k][j][id.x[1]][l];
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
  /*
  if (rank == 0) { 
    PetscPrintf(PETSC_COMM_WORLD,"step = %d\n",user->cnt);
    for (k=-1;k<=1;k++) for (j=-1;j<=1;j++) for (i=-1;i<=1;i++) 
      if (! ( (i==-1 && j==-1) || (i==-1 && k==-1) || (j==-1 && k==-1) ) ) {
        PetscPrintf(PETSC_COMM_SELF,"@[%2d][%2d][%2d] nO2+=%2.5e; nCO2+=%2.5e; nO+=%2.5e; vO2+=[%2.5e %2.5e %2.5e]; vCO2+=[%2.5e %2.5e %2.5e]; vO+=[%2.5e %2.5e %2.5e]; pO2+=%2.5e; pCO2+=%2.5e; pO+=%2.5e; pe=%2.5e; B=[%2.5e %2.5e %2.5e]\n",i,j,k,u[k][j][i][0]*user->n0,u[k][j][i][1]*user->n0,u[k][j][i][2]*user->n0,u[k][j][i][3]*user->v0,u[k][j][i][4]*user->v0,u[k][j][i][5]*user->v0,u[k][j][i][6]*user->v0,u[k][j][i][7]*user->v0,u[k][j][i][8]*user->v0,u[k][j][i][9]*user->v0,u[k][j][i][10]*user->v0,u[k][j][i][11]*user->v0,u[k][j][i][12]*user->p0,u[k][j][i][13]*user->p0,u[k][j][i][14]*user->p0,u[k][j][i][15]*user->p0,u[k][j][i][16]*user->B0,u[k][j][i][17]*user->B0,u[k][j][i][18]*user->B0);
    }
  }
  */
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

  // Get local grid boundaries //
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
  /*
  if (rank == 7) { 
    for (k=mz-2;k<=mz;k++) for (j=my-2;j<=my;j++) for (i=mx-2;i<=mx;i++) 
      if (! ( (i==mx && j==my) || (i==mx && k==mz) || (j==my && k==mz) ) ) {
        PetscPrintf(PETSC_COMM_SELF,"@[%2d][%2d][%2d] ve=[%2.5e %2.5e %2.5e]; J=[%2.5e %2.5e %2.5e]; E=[%2.5e %2.5e %2.5e]\n",i,j,k,v[k][j][i][0]*user->v0,v[k][j][i][1]*user->v0,v[k][j][i][2]*user->v0,v[k][j][i][3]*qe*user->n0*user->v0,v[k][j][i][4]*qe*user->n0*user->v0,v[k][j][i][5]*qe*user->n0*user->v0,v[k][j][i][6]*user->v0*user->B0,v[k][j][i][7]*user->v0*user->B0,v[k][j][i][8]*user->v0*user->B0);
      }
    }
    */

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

  // Get local grid boundaries //
  ierr = DMDAGetGhostCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

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
  PetscReal      /*Xmin = user->outXmin, Ymin = user->outYmin,*/ Zmin = user->outZmin;
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
  PetscReal      /*X,Y,*/Z;

  PetscFunctionBegin;

  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  // Get local grid boundaries //
  ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

  imin = 0;
  imax = mx-1;
  jmin = 0;
  jmax = my-1;
  kmin = 0;
  kmax = mz-1;

  // Calculate the flows through the boundaries of the PLOTTED domain //
  PetscPrintf(PETSC_COMM_WORLD,"START FLUX CALCULATIONS.\n");
  // W-E flows //
  for (j=ys; j<ys+ym; j++) {
    // Define y-space increment //
    if (j==0)         DY = (y[   1]-y[   0])*L    ;
    else if (j==my-1) DY = (y[my-1]-y[my-2])*L    ;
    else              DY = (y[ j+1]-y[ j-1])*L/2.0;

    for (k=zs; k<zs+zm; k++) {
      // Define z-space increment //
      if (k==0)         DZ = (z[   1]-z[   0])*L    ;
      else if (k==mz-1) DZ = (z[mz-1]-z[mz-2])*L    ;
      else              DZ = (z[ k+1]-z[ k-1])*L/2.0;

      // W-flow //
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

      // E-flow //
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

  // S-N flows //
  for (k=zs; k<zs+zm; k++) {
    // Define z-space increment //
    if (k==0)         DZ = (z[   1]-z[   0])*L    ;
    else if (k==mz-1) DZ = (z[mz-1]-z[mz-2])*L    ;
    else              DZ = (z[ k+1]-z[ k-1])*L/2.0;

    for (i=xs; i<xs+xm; i++) {

      // Define x-space increment //
      if (i==0)         DX = (x[   1]-x[   0])*L    ;
      else if (i==mx-1) DX = (x[mx-1]-x[mx-2])*L    ;
      else              DX = (x[ i+1]-x[ i-1])*L/2.0;

      // S-flow //
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
     // if(i==0) printf("@%d z=%12.6e\tnO2+=%12.6e\tvO2+=%12.6e\tDY=%12.6e\tDZ=%12.6e\n", rank, Zmin + z[k]*L, u[k][j][i][d.ni[0]]*n0,u[k][j][i][d.vi[0][1]]*v0,DY,DZ);

      // N-flow //
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

  // B-T flows //
  for (i=xs; i<xs+xm; i++) {
    // Define x-space increment //
    if (i==0)         DX = (x[   1]-x[   0])*L    ;
    else if (i==mx-1) DX = (x[mx-1]-x[mx-2])*L    ;
    else              DX = (x[ i+1]-x[ i-1])*L/2.0;

    for (j=ys; j<ys+ym; j++) {
      // Define y-space increment //
      if (j==0)         DY = (y[   1]-y[   0])*L    ;
      else if (j==my-1) DY = (y[my-1]-y[my-2])*L    ;
      else              DY = (y[ j+1]-y[ j-1])*L/2.0;

      // B-flow //
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

      // T-flow //
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

  // Summing subdomain values //
  /*
  m = MPI_Allreduce(MPI_IN_PLACE,&Fi[0][0],3*6,MPIU_REAL,MPIU_SUM,PETSC_COMM_WORLD);
  m = MPI_Allreduce(MPI_IN_PLACE,&Fe      ,6  ,MPIU_REAL,MPIU_SUM,PETSC_COMM_WORLD);
  */

  // Fill the diagnotics.dat file //
  /*
  if (rank==0) {
    sprintf(fName,"%s/diagnostics.dat",user->vName);
    nFile = fopen(fName,"w");
    if (step==0) 
      fprintf(nFile,"t           \tF.w(O2+)    \tF.e(O2+)    \tF.s(O2+)    \tF.n(O2+)    \tF.b(O2+)    \tF.t(O2+)    \tF.w(CO2+)   \tF.e(CO2+)   \tF.s(CO2+)   \tF.n(CO2+)   \tF.b(CO2+)   \tF.t(CO2+)   \tF.w(O+)     \tF.e(O+)     \tF.s(O+)     \tF.n(O+)     \tF.b(O+)     \tF.t(O+)     \tF.w(e)      \tF.e(e)      \tF.s(e)      \tF.n(e)      \tF.b(e)      \tF.t(e)      \n");
    
    fprintf(nFile,"%12.6e\t",t);
    for (l=0; l<3; l++) 
      for (m=0; m<6; m++) 
        fprintf(nFile,"%12.6e\t",Fi[l][m]);
    for (m=0; m<6; m++)
      fprintf(nFile,"%12.6e\t",Fe[m]);
    fprintf(nFile,"\n");
    fclose(nFile);
  }
  */
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
//tolerance      eps = user->eps;
  DM             da = (DM)user->da;
  DM             db = (DM)user->db;
  svi            s = user->s;
  
  PetscInt       mx = user->mx, my = user->my, mz = user->mz;
  PetscReal      n0 = user->n0, v0 = user->v0, p0 = user->p0, T0 = user->T0;
  PetscReal      un[3] = {user->un[0]/v0, user->un[1]/v0, user->un[2]/v0};       // neutral wind
  PetscInt       rank;
  PetscInt       i,j,k,m,xs,ys,zs,xm,ym,zm;
  PetscReal      vi[3][3],ve[3]; // velocity of i,e
  PetscReal      E[3];  // E-field
  PetscReal      J[3];  // total conduction current
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
  
  PetscFunctionBegin;

  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  Lz = user->outZmax - user->outZmin;

  ierr = DMDAVecGetArrayDOF(db,V,&v);CHKERRQ(ierr);
  id.x[0] = 0; id.y[0] = 0; id.z[0] = 0;
  id.x[1] = 0; id.y[1] = 0; id.z[1] = 0;
  id.x[2] = 0; id.y[2] = 0; id.z[2] = 0;

  dh.x[0] = 0; dh.y[0] = 0; dh.z[0] = 0;
  dh.x[1] = 0; dh.y[1] = 0; dh.z[1] = 0;

  // Get local grid boundaries //
  ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

//printf("@%d, xs=%d\txs+xm=%d\tys=%d\tys+ym=%d\tzs=%d\tzs+zm=%d\n",rank,xs,xs+xm,ys,ys+ym,zs,zs+zm);
//exit(1);

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
      Y = Ymin + y[j]*L;
      for (i=xs; i<xs+xm; i++) {
        X = Xmin + x[i]*L;
        if (bcType == 0 || bcType == 1 || bcType == 2 || bcType == 10 || bcType == 11) {
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


        /*
        Te          = pe*p0/(ne*n0*kB);
        nn[CO2]     = Interpolate(user->RefProf, 5 ,Z, lin_flat);
        nn[O]       = Interpolate(user->RefProf, 4 ,Z, lin_flat);
        nu[CO2]     = Ven(CO2, nn[CO2], Te);
        nu[O]       = Ven(O  , nn[O]  , Te);
        */

        Te          = pe/ne;
        nn[CO2]     = Interpolate(user->RefProf, 5 ,Z, lin_flat)/n0;
        nn[O]       = Interpolate(user->RefProf, 4 ,Z, lin_flat)/n0;
        nu[0]       = v41(nn[CO2] *n0 , Te      *T0);                // e    CO2
        nu[1]       = v42(nn[O]   *n0 , Te      *T0);                // e    O

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

        // Current //
        J[0] = dB[2].dy - dB[1].dz; 
        J[1] = dB[0].dz - dB[2].dx;
        J[2] = dB[1].dx - dB[0].dy;
       
        // Electron velocity //
        ve[0] = ( (ni[O2p]*vi[O2p][0] + ni[CO2p]*vi[CO2p][0] + ni[Op]*vi[Op][0]) - J[0] )/ ne;
        ve[1] = ( (ni[O2p]*vi[O2p][1] + ni[CO2p]*vi[CO2p][1] + ni[Op]*vi[Op][1]) - J[1] )/ ne;
        ve[2] = ( (ni[O2p]*vi[O2p][2] + ni[CO2p]*vi[CO2p][2] + ni[Op]*vi[Op][2]) - J[2] )/ ne;
        /*
        for (m=0; m<3; m++)
          if (ne<=user->eps.n){
            ve[m] = 0.0;
            PetscPrintf(PETSC_COMM_WORLD,"HERE is the PROBLEM\n");
          }
        */
        //if (i==0 && j==0 && k==0) printf("dz=%2.1e km, nu=%2.5e ve=[%2.3e %2.3e %2.3e] un=[%2.3e %2.3e %2.3e]\n", dh.z[0],nu[CO2]+nu[O], ve[0]*v0,ve[1]*v0,ve[2]*v0,un[0]*v0,un[1]*v0,un[2]*v0);
        
        // E-field //
        if (vDamping) { //((1.0-tanh(Z-Lz)/(Lz/12.0))/2.0);
          for (m=0; m<3; m++) { 
            un[m]*=
              .5*(1+erf((X-XL)/lambda)) * .5*(1-erf((X-XU)/lambda)) *
              .5*(1+erf((Y-YL)/lambda)) * .5*(1-erf((Y-YU)/lambda)) *
              .5*(1+erf((Z-ZL)/lambda)) * .5*(1-erf((Z-ZU)/lambda)) ;
          }
        }
        E[0] = -CrossP(ve,B,0) - dpe.dx/ne + (nu[CO2] + nu[O]) * tau * (un[0] - ve[0]);
        E[1] = -CrossP(ve,B,1) - dpe.dy/ne + (nu[CO2] + nu[O]) * tau * (un[1] - ve[1]);
        E[2] = -CrossP(ve,B,2) - dpe.dz/ne + (nu[CO2] + nu[O]) * tau * (un[2] - ve[2]);

        // Store data //
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
                 };                                                              // * at the CURRENT point (k,j,i)
  
  PetscReal      vA_max[3] = {                                                   // Maximum Alfven speed
                    0.0,                                                         // * 
                    0.0,                                                         // * for a GIVEN species l
                    0.0                                                          // * in ANY direction m
                 };                                                              // * at the ENTIRE SUBDOMAIN
  
  PetscReal      vs[3] = {                                                       // Sonic speed
                    0.0,                                                         // * 
                    0.0,                                                         // * for a GIVEN species l
                    0.0                                                          // * in ANY direction m
                 };                                                              // * at the CURRENT point (k,j,i)

  PetscReal      vs_max[3] = {                                                   // Maximum Sonic speed
                    0.0,                                                         // * 
                    0.0,                                                         // * for a GIVEN species l
                    0.0                                                          // * in ANY direction m
                 };                                                              // * at the ENTIRE SUBDOMAIN
  
  PetscReal      vf[3] = {                                                       // Fast magnetosonic wave speed
                    0.0,                                                         // * 
                    0.0,                                                         // * for a GIVEN species l
                    0.0                                                          // * in ANY direction m
                 };                                                              // * at the CURRENT point (k,j,i)

  PetscReal      vf_max[3] = {                                                   // Maximum Fast magnetosonic wave speed
                    0.0,                                                         // * 
                    0.0,                                                         // * for a GIVEN species l
                    0.0                                                          // * in ANY direction m
                 };                                                              // * at the ENTIRE SUBDOMAIN

  PetscReal      ve_max[3] = {                                                   // Maximum speed of the hydrodynamic electron flow
                    0.0,                                                         // * 
                    0.0,                                                         // * 
                    0.0                                                          // * in ANY direction m
                 };                                                              // * at the ENTIRE SUBDOMAIN

  PetscReal      vcr[3][3] = {                                                   // Critical speed
                    {0.0, 0.0, 0.0},                                             // * sum of the hydrodynamic and fast magnetoacoustic wave speeds
                    {0.0, 0.0, 0.0},                                             // * for a GIVEN species l
                    {0.0, 0.0, 0.0}                                              // * in the GIVEN direction m
                 };                                                              // * at the CURRENT point (k,j,i)

  PetscReal      vcr_max[3][3] = {                                               // Critical speed 
                    {0.0, 0.0, 0.0},                                             // * sum of the hydro-dynamic and fast magnetoacoustic wave speeds
                    {0.0, 0.0, 0.0},                                             // * for a GIVEN species l
                    {0.0, 0.0, 0.0}                                              // * in the GIVEN direction m
                 };                                                              // * in the SUBDOMAIN (at the end of the space scanning of the subdomain)

  PetscReal      vth = 0.0;                                                      // Thermal velocity

  PetscReal      T_cr[4] = {0.0,0.0,0.0,0.0};
  PetscErrorCode ierr;
  PetscReal      t,dt,dt_CFL,h[3];
  PetscReal      nu[4][6];                                                       // ion-neutral collision frequency
  PetscReal      E[3],B[3];                                                      // E-field, B-field
  PetscReal      ne,ni[3],nn[2];                                                 // [O2+], [CO2+], [O+], ne in m^-3
  PetscReal      Ti[3],Te,Tn[2];                                                 // Temperature of e, O2+, CO2+, O+, CO2, 0 in K
  PetscReal      pe,pi[3];                                                       // pressure of O2+, CO2+, O+, e in J/m^3
  Vec            V;
  PetscReal      ****u,****v,****f;
  Vec            localU,localV;
  PetscInt       rank,step;
  PetscInt       i,j,k,l,m,n,xs,ys,zs,xm,ym,zm;
  PetscReal      ve[3],vi[3][3];                                                 // electron and ion velocities 
  diff           dE[3],dpe,dve[3],dni[3],dpi[3],dvi[3][3];                       // Other derivatives 
  stencil        id;
  coeff3         D1;
  coeff2         dh;
  PetscBool      flag = PETSC_FALSE;
  
  PetscFunctionBegin;
  id.x[0] = 0; id.y[0] = 0; id.z[0] = 0;
  id.x[1] = 0; id.y[1] = 0; id.z[1] = 0;
  id.x[2] = 0; id.y[2] = 0; id.z[2] = 0;

  dh.x[0] = 0; dh.y[0] = 0; dh.z[0] = 0;
  dh.x[1] = 0; dh.y[1] = 0; dh.z[1] = 0;
  
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  Lz = user->outZmax - user->outZmin;

  // Get current iteration t, dt, step //
  ierr = TSGetTimeStep(ts,&dt);CHKERRQ(ierr);
  ierr = TSGetTime(ts,&t);CHKERRQ(ierr);
  ierr = TSGetTimeStepNumber(ts,&step);CHKERRQ(ierr);

  ierr = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr = DMGetGlobalVector(db,&V);CHKERRQ(ierr);
  
  // Create and fill ****u inner cells proc ghosts //
  ierr = DMGetLocalVector(da,&localU);CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(da,U,INSERT_VALUES,localU);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da,U,INSERT_VALUES,localU);CHKERRQ(ierr);

  // Fill ****u domain ghosts //
  ierr = DMDAVecGetArrayDOF(da,localU,&u);CHKERRQ(ierr);
  ierr = FormBCu(u,user);CHKERRQ(ierr);

  // Create and fill ****v inner cells//
  ierr = FormIntermediateFunction(u,V,user);CHKERRQ(ierr);

  // Create and fill ****v inner cells proc ghosts //
  ierr = DMGetLocalVector(db,&localV);CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(db,V,INSERT_VALUES,localV);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(db,V,INSERT_VALUES,localV);CHKERRQ(ierr);
  
  // Fill ****v domain ghosts //
  ierr = DMDAVecGetArrayDOF(db,localV,&v);CHKERRQ(ierr);
  ierr = FormBCv(v,user);CHKERRQ(ierr);

  // Check values sanity //
  ierr = CheckValues(u,v,user);CHKERRQ(ierr);

  // Form F(u,v) = Udot //
  ierr = DMDAVecGetArrayDOF(da,F,&f);CHKERRQ(ierr);

  // Compute function over the locally owned part of the grid //
  dt_CFL = 1e300;
  flag   = PETSC_FALSE;
  for (k=zs; k<zs+zm; k++) {
    Z = Zmin + z[k]*L;
    for (j=ys; j<ys+ym; j++) {
      Y = Ymin + y[j]*L;
      for (i=xs; i<xs+xm; i++) {
        X = Xmin + x[i]*L;
        if (bcType == 0 || bcType == 1 || bcType == 2 || bcType == 10 || bcType == 11) {
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

        /*
        nn[CO2] = Interpolate(user->RefProf, 5,Z,lin_exp ); 
        nn[O]   = Interpolate(user->RefProf, 4,Z,lin_exp );
        Tn[CO2] = Interpolate(user->RefProf,10,Z,lin_flat);
        Tn[O]   = Interpolate(user->RefProf,10,Z,lin_flat);
        for (l=0; l<3; l++) {
          Ti[l] = pi[l]*p0/(ni[l]*n0*kB);
          for (n=0; n<2; n++) {
            nu[l][n] = Vin(l,n,nn[n],Ti[l],Tn[n]);
          }
        }
        Te = pe*p0/(ne*n0*kB);
        for (n=0; n<2; n++) 
          nu[e][n] = Ven(n,nn[n],Te);

        for (l=0; l<4; l++)
          for (n=0; n<2; n++)
            PetscPrintf(PETSC_COMM_WORLD,"nu[%d,%d]=%e\n",l,n,nu[l][n]);
        */
        nn[CO2] = Interpolate(user->RefProf, 5,Z,lin_exp )/n0; 
        nn[O]   = Interpolate(user->RefProf, 4,Z,lin_exp )/n0;
        for (l=0; l<3; l++)
          Ti[l] = pi[l]/ni[l];
        Tn[CO2] = Interpolate(user->RefProf,10,Z,lin_flat)/T0;
        Tn[O]   = Interpolate(user->RefProf,10,Z,lin_flat)/T0;
        Te      = pe/ne;

        // collisions //
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

        /*
        for (l=0; l<4; l++)
          for (n=0; n<2; n++)
            PetscPrintf(PETSC_COMM_WORLD,"nu[%d,%d]=%e\n",l,n,nu[l][n]);
        exit(1);
        */
         
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

        // CFL //
        h[0] = MinAbs(dh.x[0],dh.x[1]);
        h[1] = MinAbs(dh.y[0],dh.y[1]);
        h[2] = MinAbs(dh.z[0],dh.z[1]);

        // Critical speed and time for ions //
        for (l=0; l<3; l++) {
          vs[l] = PetscSqrtScalar( (gama[e]*pe + gama[l]*pi[l]) / (ni[l]*mi[l]/me) );
          vs_max[l] = MaxAbs( vs_max[l], vs[l] );

          vA[l] = PetscSqrtScalar( Norm2(B) / (ni[l]*mi[l]/me) );
          vA_max[l] = MaxAbs( vA_max[l], vA[l] );

          vf[l] = c*PetscSqrtScalar( (vs[l]*vs[l]+vA[l]*vA[l]) / (c*c+vA[l]*vA[l]) );
          vf_max[l] = MaxAbs( vf_max[l], vf[l] );

          for (m=0; m<3; m++) {
            vcr[l][m] = PetscAbsScalar(vi[l][m]) + vf[l];
            vcr_max[l][m] = MaxAbs( vcr_max[l][m], vcr[l][m] );
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
        // Critical speed and time for electrons //
        for (m=0; m<3; m++) {
          ve_max[m] = MaxAbs(MaxAbs(ve_max[m], ve[m]),vth);
          for (l=0; l<3; l++) {
            if (ve_max[m] > vcr_max[l][m]) { 
              /* 
               * Case when the electron related velocities are the critical velocities
               */
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

        /*
        if(vp > vA) PetscPrintf(PETSC_COMM_SELF, "SuperAlfvenic speed @ [%d,%d,%d], |vloc|=%2.3e, vA=%2.3e (Z=%f km)\n",i,j,k,vp*v0,vA*v0, Z*1e-3);
        assert(vp <= v,,A); // Check sub-Alfvenic speeds
        if(vp > vs_n) PetscPrintf(PETSC_COMM_SELF, "Supersonic (neutral) speed @ [%d,%d,%d], |vloc|=%2.3e, vs_n=%2.3e (Z=%f km)\n",i,j,k,vp*v0,vs_n*v0, Z*1e-3);
        assert(vp <= vs_n); // Check subsonic speeds
        if(vp > vs_c) PetscPrintf(PETSC_COMM_SELF, "Supersonic (charged) speed @ [%d,%d,%d], |vloc|=%2.3e, vs_c=%2.3e (Z=%f km)\n",i,j,k,vp*v0,vs_c*v0, Z*1e-3);
        assert(vp <= vs_c); // Check subsonic speeds
        */

        // Equations //
        // Eq 0-2: Continuity equations for ions //
        for (l=0; l<3; l++)
          f[k][j][i][d.ni[l]] = -ni[l]*(dvi[l][0].dx + dvi[l][1].dy + dvi[l][2].dz) - (vi[l][0]*dni[l].dx + vi[l][1]*dni[l].dy + vi[l][2]*dni[l].dz);
      
        // Eq 3-11: Ion momenta // NEGLECT O and e-n collisions FOR NOW (07-29-11) //
        if (vDamping) { //((1.0-tanh(Z-Lz)/(Lz/12.0))/2.0);
          for (m=0; m<3; m++) { 
            un[m]*=
              .5*(1+erf((X-XL)/lambda)) * .5*(1-erf((X-XU)/lambda)) *
              .5*(1+erf((Y-YL)/lambda)) * .5*(1-erf((Y-YU)/lambda)) *
              .5*(1+erf((Z-ZL)/lambda)) * .5*(1-erf((Z-ZU)/lambda)) ;
          }
        }

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
      }
    }
  }

  m = MPI_Allreduce(MPI_IN_PLACE,&vs_max[0]  ,3 ,MPIU_REAL,MPIU_MAX,PETSC_COMM_WORLD);
  m = MPI_Allreduce(MPI_IN_PLACE,&vA_max[0]  ,3 ,MPIU_REAL,MPIU_MAX,PETSC_COMM_WORLD);
  m = MPI_Allreduce(MPI_IN_PLACE,&vf_max[0]  ,3 ,MPIU_REAL,MPIU_MAX,PETSC_COMM_WORLD);
  m = MPI_Allreduce(MPI_IN_PLACE,&ve_max[0]  ,3 ,MPIU_REAL,MPIU_MAX,PETSC_COMM_WORLD);
  m = MPI_Allreduce(MPI_IN_PLACE,&vcr_max[0] ,9 ,MPIU_REAL,MPIU_MAX,PETSC_COMM_WORLD);

  /*
  for (m=0; m<3; m++) {
    for (l=0; l<3; l++) {
      if (ve_max[m] > vcr_max[l][m]) {
        if(blockers) {
          PetscPrintf(PETSC_COMM_WORLD,"The CFL criterion is not stringent enough!\n");
          exit(12);
        }
      }
    }
  }
  */

  /* dt_CFL:
   * critical time for CFL condition
   * in the ENTIRE domain (i.e., for all point (k,j,i))
   * for ALL species (usually denoted l)
   * indepently from ANY direction in space (usually denoted m)
   */
  m = MPI_Allreduce(MPI_IN_PLACE,&dt_CFL ,1 ,MPIU_REAL,MPIU_MIN,PETSC_COMM_WORLD);

  // Adapt timestep to satisfy CFL //
  user->dt = c_CFL*dt_CFL;

  /*
   * BELOW IS THE LIST OF DIAGNOSTIC TO ACTIVATE
  ierr = PetscPrintf(PETSC_COMM_WORLD,
    "\t==> t=%12.5e\t new dt=%12.5e\n\t==> O2+  max velocities: vs=%+14.7e, vA=%+14.7e, vf=%+14.7e, vcr=[%+14.7e,%+14.7e,%+14.7e]\n\t==> CO2+ max velocities: vs=%+14.7e, vA=%+14.7e, vf=%+14.7e, vcr=[%+14.7e,%+14.7e,%+14.7e]\n\t==> O+   max velocities: vs=%+14.7e, vA=%+14.7e, vf=%+14.7e, vcr=[%+14.7e,%+14.7e,%+14.7e]\n",
     ftime*tau,user->dt*tau,
     vs_max[O2p] *v0, vA_max[O2p] *v0, vf_max[O2p] *v0, vcr_max[O2p][0] *v0 ,vcr_max[O2p][1] *v0 ,vcr_max[O2p][2] *v0,
     vs_max[CO2p]*v0, vA_max[CO2p]*v0, vf_max[CO2p]*v0, vcr_max[CO2p][0]*v0 ,vcr_max[CO2p][1]*v0 ,vcr_max[CO2p][2]*v0,
     vs_max[Op]  *v0, vA_max[Op]  *v0, vf_max[Op]  *v0, vcr_max[Op][0]  *v0 ,vcr_max[Op][1]  *v0 ,vcr_max[Op][2]  *v0);
   */
  // Calculates and write fluxes at the domain boundaries //
  /*
  if(step%user->viz_dstep==0) {
    if(user->fluxRec) {
      ierr = CalculateFluxes(t,step,u,v,user);CHKERRQ(ierr);
      user->fluxRec = PETSC_FALSE;
    } else {
      user->fluxRec = PETSC_TRUE;
    }
    PetscPrintf(PETSC_COMM_WORLD,"%.6i %i\n",step,user->fluxRec);
  }
  */


  // Restore vectors //
  ierr = DMDAVecRestoreArrayDOF(db,localV,&v);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(da,localU,&u);CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(da,F,&f);CHKERRQ(ierr);
  
  ierr = DMRestoreLocalVector(db,&localV);CHKERRQ(ierr);
  //ierr = DMLocalToGlobalBegin(da,localU,INSERT_VALUES,U);CHKERRQ(ierr);
  //ierr = DMLocalToGlobalEnd(da,localU,INSERT_VALUES,U);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(da,&localU);CHKERRQ(ierr);

  // Destroy V //
  ierr = DMRestoreGlobalVector(db,&V); CHKERRQ(ierr);
  PetscFunctionReturn(0); 
} 

