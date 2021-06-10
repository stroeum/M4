/*
 * MyCtx.h
 * Created by Jeremy A Riousset on 02/04/11
 * User defined functions and structures for use with DAE solvers
 */

#define eps0 8.85418782e-12    //_F/_m, free space permittivity
#define mu0 1.2566370614e-6    //_N/_A, free space permeability
#define kB 1.3806503e-23       //_J/_K, Boltzmann constant
#define qe 1.60217646e-19      //_C, elementary charge
#define G 6.673e-11            //_m^3/_kg/_s^/\2, gravitational constant
#define c_CFL 0.5              //_, CFL coefficient 
#define Nz_REF 301             // number of reference altitude points
#define Np_REF 18              // number of parameters in Profiles.dat
#define MAX_LINE_LENGTH 1024   // maximum number of characters in a line

#ifndef MYCTX_H
#define MYCTX_H

#include <petsctime.h>
#include <petscdmda.h>
#include <petscts.h>
#include <petsc/private/tsimpl.h>
#include <unistd.h>
#include <assert.h>
#include "MyProfiles.h"

enum chargedspecies {O2p, CO2p, Op, e};
enum neutrals {CO2, O};

typedef struct {
	PetscErrorCode (*onestep)(TS,PetscReal,PetscReal,Vec);
	char *type_name;
	PetscInt nstages;
	Vec *work;
	PetscInt nwork;
	PetscBool workout;
} TS_SSP;

typedef struct{
	PetscReal n,p,v;
} tolerance;

typedef struct{
	PetscReal dx,dy,dz;
} diff;

typedef struct{
	PetscInt x[3],y[3],z[3];
} stencil;

typedef struct{
	PetscReal x[3],y[3],z[3];
} coeff3;

typedef struct{
	PetscReal x[2],y[2],z[2];
} coeff2;

typedef struct{
	PetscInt ni[3],vi[3][3],pi[3],pe,B[3];
} dvi;

typedef struct{
	PetscInt ve[3],J[3],E[3];
} svi;

typedef struct{
	PetscReal g[3],dp[4][3],E[4][3],B[4][3],col[4][3],adv[4][3];
} fdv; // variables with force dependent terms

/* MonitorCtx: used by MyTSMonitor() */
typedef struct {
	PetscBool drawcontours;   
} MonitorCtx;

/* AppCtx: used by FormFunction() and FormJacobian() */
typedef struct {
	PetscLogDouble t0;                            // Simulation start _s
	tolerance      eps;                           // Minimum value for densities
	PetscInt       cnt;                           // A counter
	PetscInt       mx,my,mz;                      // Nb of gridpoints in the x, y, and z-directions for the simulation domain
	PetscReal      *x,*y,*z;                      // x-, y-, and z-vectors
	PetscReal      inXmin,inYmin,inZmin;          // x-, y-, and z-minima
	PetscReal      inXmax,inYmax,inZmax;          // x-, y-, and z-maxima
	PetscReal      outXmin,outYmin,outZmin;       // x-, y-, and z-minima
	PetscReal      outXmax,outYmax,outZmax;       // x-, y-, and z-maxima
	PetscReal      vizbox[6];                     // Visualization box
	PetscInt       viz_dstep;                     // Steps skipped for vizualization
	PetscReal      rM,mM,gM;                      // Mars' radius, mass, and gravitational acceleration
	PetscReal      me,mi[3];                      // me,mi[O2+,CO2+,O+]
	PetscReal      Lx,Ly,Lz,L;                    // Domain size in x-,y-, and z-directions, and reference length
	PetscReal      dt;                            // time-step
	PetscReal      ti,tf;                         // initial time and final time
	PetscInt       istep,maxsteps;                // initial step, maximum number of time iterations
	PetscReal      tau;                           // L/vo
	PetscBool      xtra_out;                      // output all variables
	PetscBool      smoothing;                     // use smoothing Y/N
	PetscBool      blockers;                      // use blockers Y/N
	PetscBool      limiters;                      // use smoothing Y/N
	PetscBool      vDamping;                      // use damping to get null velocities at 400 km
	PetscReal      ZL,ZU,lambda;                  // damping parameters
	PetscInt       BfieldType;                    // Define the B-field configuration a t=0 
	PetscReal      B[4];                          // x-, y-, and z-component of the default uniform magnetic field
	PetscReal      un[3];                         // x-, y-, and z-component of the neutral wind
	PetscReal      ui[3];                         // Initial x-, y-, and z-component of the ions
	PetscReal      n0;                            // no[O2+,CO2+,O+,e], reference particule number densities
	PetscReal      gama[4];                       // Specific heat ratio (O2+,CO2+,O+,e)
	PetscReal      v0;                            // Reference velocity (thermal velocity)
	PetscReal      p0;                            // Reference kinectic pressure
	PetscReal      T0;                            // Reference temperature
	PetscReal      B0;                            // Reference magnetic field
	svi            s;                             // Definition of indices for "static" variables (defined by an algebraic equation)
	dvi            d;                             // Definition of indices for "dynamic" variables (defined by a partial differential equation)
	PetscBool      isInputFile;                   // Input file: 1, and 0, no input file
	PetscBool      bcRec;                         // Record BC
	char           InputFile[PETSC_MAX_PATH_LEN]; // Name of the input file
	char           dName[PETSC_MAX_PATH_LEN];     // Name of the output directory
	char           vName[PETSC_MAX_PATH_LEN];     // Name of the visualization direction
	PetscInt       ldName;
	DM             da,db;
	PetscInt       jacType;
	PetscInt       bcType;
	PetscInt       chemswitch;                    // Turn chemistry on/off
	PetscInt       collswitch;                    // Turn collisions on/off
	PetscInt       gravswitch;                    // Turn gravity on/off
	PetscReal      RefProf[Np_REF][Nz_REF];       // Table of reference Profiles
	PetscReal      RefPart[4][Nz_REF];            // Table of reference Partition
} AppCtx;

extern PetscErrorCode InitCtx(AppCtx*,MonitorCtx*);

extern PetscErrorCode OutputData(void*);

extern PetscErrorCode CFL(TS);

extern PetscReal      CrossP(PetscReal*,PetscReal*,PetscInt);
extern PetscReal      SumAbs(PetscReal,PetscReal);
extern PetscReal      MaxAbs(PetscReal,PetscReal);
extern PetscReal      MinAbs(PetscReal,PetscReal);
extern PetscReal      V_Dipole(PetscReal mu, PetscReal xs, PetscReal ys, PetscReal zs, PetscReal x, PetscReal y, PetscReal z, PetscInt m);
extern PetscReal      H_Dipole(PetscReal mu, PetscReal xs, PetscReal ys, PetscReal zs, PetscReal x, PetscReal y, PetscReal z, PetscInt m);
extern PetscReal      Arcades(PetscReal x, PetscReal y, PetscReal z, PetscInt m);
extern PetscReal      MultiArcades(PetscReal x, PetscReal y, PetscReal z, PetscInt m);
extern PetscReal      Norm2(PetscReal*); // Norm 2 squared!

extern PetscErrorCode TSSSPStep_RK_2_JAR(TS,PetscReal,PetscReal,Vec);
extern PetscErrorCode TSSSPStep_LAX(TS,PetscReal,PetscReal,Vec);
extern PetscErrorCode TSSSPStep_LW(TS,PetscReal,PetscReal,Vec);
extern PetscErrorCode WeightedAverage(Vec,void*);
extern PetscErrorCode Average(Vec,void*);

#endif

