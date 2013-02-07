/*
 * MyProfiles.h
 * Created by Jeremy A Riousset & Carol S Paty on 03/16/11
 *
 * Define a set of very useful profiles between 50 and 450 km.
 */

#include "MyProfiles.h"
#define alpha 3.0

#undef __FUNCT__
#define __FUNCT__ "ReadTable"
PetscInt  ReadTable(PetscReal table[][Nz_REF], PetscInt P, const char *fName)
{
  char buffer[MAX_LINE_LENGTH];
  FILE *pFile;
  int m=0, p=0;

  pFile=fopen(fName,"r");
  fgets(buffer, MAX_LINE_LENGTH, pFile);
  
  for (m=0; m<Nz_REF; m++) {
    for (p=0; p<P; p++) {
      fscanf(pFile,"%lf", &table[p][m]);
    }
  }
  fclose(pFile);
    
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "Vin"

/* ------------------------------------------------------------------- */
/* 
   nu = Rate(i,n,ni,nn,Ti,Tn,h) - return the collision frequency vin at the altitude h

   Input Parameters:
.  i - choice of the (c)harged species
      0. O2+
      1. CO2+
      2. O+
.  n - choice of (n)eutral species
      0. CO2
      1. O
.  nc - charged species density (m^-3)
.  nn - neutral density (m^-3)
.  Tn - neutral gas temperature (K)
.  Tc - charged species temperature (K)
.  h - altitude
   Output Parameter:
.  nu - output value 
*/
PetscReal Vin(PetscInt i, PetscInt n, PetscReal nn, PetscReal Ti, PetscReal Tn)
{
  PetscReal  Tr,nu;

  /* Formulas from Schunk and Nagy, 2000, pp. 97--99, Tables 4.4, 4.5, 4.6 */
  Tr = (Ti+Tn)/2.0;
  //if (c==CO2p && n==CO2) printf( "Tr = %f\n", Tr);
  //if (c==Op   && n==O)   printf( "Tr = %f\n", Tr);
  if      (i==O2p  && n==CO2) {nu = 5.63e-10*nn*1e-6;}
  else if (i==CO2p && n==CO2) {nu = (Tr>850)*2.85e-11*(nn*1e-6)*sqrt(Tr)*(1-0.083*log10(Tr))*(1-0.083*log10(Tr));}
  else if (i==Op   && n==CO2) {nu = 8.95e-10*nn*1e-6;}
  else if (i==O2p  && n==O)   {nu = 2.31e-10*nn*1e-6;}
  else if (i==CO2p && n==O)   {nu = 1.76e-10*nn*1e-6;}
  else if (i==Op   && n==O)   {nu = (Tr>235)*3.67e-11*(nn*1e-6)*sqrt(Tr)*(1-0.064*log10(Tr))*(1-0.064*log10(Tr));}
  else    {printf("Warning: unknown collision\n\tAssuming null collision frequency.\n"); nu=0.0;}

  if(!(nu>=0)) nu = 0.0; //printf("h=%f, n=%d,Te=%f,Tn=%f, nn=%f,ne=%f,nu=%f\n",h,n,Te,Tn,nn,ne,nu);

  return nu;
}

#undef __FUNCT__
#define __FUNCT__ "vin"

/* ------------------------------------------------------------------- */
/* 
   nu = vin(c,n,h) - return the collision frequency vin at the altitude h

   Input Parameters:
.  c - choice of the (c)harged species
      0. O2+
      1. CO2+
      2. O+
.  n - choice of (n)eutral species
      0. CO2
      1. 0
.  h - altitude
   Output Parameter:
.  nu - output value 
*/
PetscReal vin(PetscReal profiles[][Nz_REF],PetscInt i, PetscInt n, PetscReal h)
{
  PetscReal nn=0.0,Ti,Tn;

  if (i==O2p || i==CO2p || i==Op) {
    Ti = Interpolate(profiles,9,h,lin_flat);
  } else {
    printf("Error: Unknown charged species\n");
    exit(0);
  }

  if (n==CO2) nn = Interpolate(profiles,5,h,lin_exp);
  if (n==O)   nn = Interpolate(profiles,4,h,lin_exp);
  Tn = Interpolate(profiles,10,h,lin_flat);

  return Vin(i,n,nn,Ti,Tn);
}

#undef __FUNCT__
#define __FUNCT__ "Ven"

/* ------------------------------------------------------------------- */
/* 
   nu = Rate(n,ne,nn,Te,Tn,h) - return the collision frequency ven at the altitude h

   Input Parameters:
.  n - choice of (n)eutral species
      0. CO2
      1. 0
.  ne - electron density (m^-3)
.  nn - neutral density (m^-3)
.  Te - electron temperature (K)
.  Tn - neutral gas temperature (K)
.  h - altitude
   Output Parameter:
.  nu - output value 
*/
PetscReal Ven(PetscInt n, PetscReal nn, PetscReal Te)
{
  PetscReal  nu;

  /* Formulas from Schunk and Nagy, 2000, pp. 97--99, Tables 4.4, 4.5, 4.6 */
  if      (n==CO2) {nu = 3.68e-8*(nn*1e-6)*(1+4.1e-11*pow(fabs(4500 - Te),2.93));}
  else if (n==O)   {nu = 8.9e-11*(nn*1e-6)*(1+5.7e-4*Te)*sqrt(Te);}
  else    {printf("Warning: unknown collision\n\tAssuming null collision frequency.\n"); nu=0.0;}

  if(!(nu>=0)) nu = 0.0; //printf("h=%f, n=%d,Te=%f,Tn=%f, nn=%f,ne=%f,nu=%f\n",h,n,Te,Tn,nn,ne,nu);

  return nu;
}

#undef __FUNCT__
#define __FUNCT__ "ven"

/* ------------------------------------------------------------------- */
/* 
   nu = ven(n,h) - return the collision frequency ven at the altitude h

   Input Parameters:
.  n - choice of (n)eutral species
      0. CO2
      1. 0
.  h - altitude
   Output Parameter:
.  nu - output value 
*/
PetscReal ven(PetscReal profiles[][Nz_REF],PetscInt n, PetscReal h)
{
  PetscReal nn=0.0,Te;

  Te = Interpolate(profiles,8,h,lin_flat);

  if (n==CO2) nn = Interpolate(profiles,5,h,lin_exp);
  if (n==O)   nn = Interpolate(profiles,4,h,lin_exp);

  return Ven(n,nn,Te);
}

#undef __FUNCT__
#define __FUNCT__ "Interpolate"
PetscReal Interpolate(PetscReal profiles[][Nz_REF], PetscInt i, PetscReal h, PetscInt ItpType)
{
  PetscInt   id[3];
  PetscReal  v = 0.0, Du = 0.0, Dz[3] = {0.0, 0.0, 0.0}, dh[2] = {0.0,0.0},  H = 0.0;
  PetscReal  z[3],u[3];
 
  /*
  if(h<100 || h>400) {
    printf("h is outside of the 100--400 km interval\n");
    exit(0);
  }*/

  if(h>=profiles[0][Nz_REF-1] && h<=profiles[0][0]) {
    id[2] = 0;
    id[0] = Nz_REF-1;
    while (id[2]!=id[0]-1)
    {
      id[1] = (PetscInt)(id[0]+id[2])/2;
      dh[0] = h-profiles[0][id[1]];
      if(dh[0]<0)  id[2]=id[1];
      if(dh[0]>=0) id[0]=id[1];
    };

    z[0] = profiles[ 0 ][id[0]];
    z[1] = profiles[ 0 ][id[2]];
    u[0] = profiles[i+1][id[0]];
    u[1] = profiles[i+1][id[2]];
    /*
    if(u[0]!=0)
      if(u[1]/u[0]>0) // Exponential interpolation
        v = u[0]*pow(u[1]/u[0],(h-z[0])/(z[1]-z[0]));
      else  // Linear interpolation if change of sign
        v = (h-z[0])/(z[1]-z[0])*u[1] + (h-z[1])/(z[0]-z[1])*u[0];
    else
      v = 0;
    */
    v = (h-z[0])/(z[1]-z[0])*u[1] + (h-z[1])/(z[0]-z[1])*u[0];
  }
  else if (h>profiles[0][0]) {
    id[0] = 0;
    id[1] = 1;
    id[2] = 2;
    z[0]  = profiles[ 0 ][id[0]];
    z[1]  = profiles[ 0 ][id[1]];
    z[2]  = profiles[ 0 ][id[2]];
    u[0]  = profiles[i+1][id[0]];
    u[1]  = profiles[i+1][id[1]];
    u[2]  = profiles[i+1][id[2]];
    dh[0] = z[1] - z[0];
    dh[1] = z[2] - z[0];
    Dz[0] = - 1.0/dh[0] - 1.0/dh[1]; 
    Dz[1] =   1.0/dh[0] - 1.0/(dh[0]-dh[1]); 
    Dz[2] =   1.0/dh[1] + 1.0/(dh[0]-dh[1]);

    Du    = Dz[0]*u[0] + Dz[1]*u[1] + Dz[2]*u[2];

    // linear interpolation above max reference altitude //
    if (ItpType==lin_lin) {
      v = Du*(h-z[0]) + u[0];
    }
    // exponential decrease above max reference altitude //
    if (ItpType==lin_exp) {
      H     = u[0]/Du;
      if (H<0) {
        if (fabs((h-z[0])/H)<alpha)
          v = u[0]*exp( (h-z[0])/H );
        else
          v = u[0]*exp(-alpha);
      } else {
        PetscPrintf(PETSC_COMM_WORLD,"Wrong Interpolation\n");
        exit(1);
      }
    }
    // flat interpolation above max reference altitude //
    if (ItpType==lin_flat) {
      v = u[0];
    }
  }
  else if (h<profiles[0][Nz_REF-1]) {
    id[0] = Nz_REF-1;
    id[1] = Nz_REF-2;
    id[2] = Nz_REF-3;
    z[0]  = profiles[ 0 ][id[0]];
    z[1]  = profiles[ 0 ][id[1]];
    z[2]  = profiles[ 0 ][id[2]];
    u[0]  = profiles[i+1][id[0]];
    u[1]  = profiles[i+1][id[1]];
    u[2]  = profiles[i+1][id[2]];
    dh[0] = z[1] - z[0];
    dh[1] = z[2] - z[0];
    Dz[0] = - 1.0/dh[0] - 1.0/dh[1]; 
    Dz[1] =   1.0/dh[0] - 1.0/(dh[0]-dh[1]); 
    Dz[2] =   1.0/dh[1] + 1.0/(dh[0]-dh[1]);

    Du    = Dz[0]*u[0] + Dz[1]*u[1] + Dz[2]*u[2];

    // linear interpolation below min reference altitude //
    if (ItpType==lin_lin) {
      v = Du*(h-z[0]) + u[0];
    }
    // exponential decrease below min reference altitude //
    if (ItpType==lin_exp) {
      H     = u[0]/Du;
      if (H<0) {
        if (fabs((h-z[0])/H)<alpha)
          v = u[0]*exp( (h-z[0])/H );
        else
          v = u[0]*exp(-alpha);
      } else {
        PetscPrintf(PETSC_COMM_WORLD,"Wrong Interpolation\n");
        exit(1);
      }
    }
    // flat interpolation below min reference altitude //
    if (ItpType==lin_flat) {
      v = u[0];
    }
  }

  return v;
}
