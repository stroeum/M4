/*
 * MyCollisions.c
 * Created by Jeremy A Riousset & Carol S Paty on 03/16/11
 *
 * Define collision frequencies between ion, electron, and neutral particles in the Martian ionosphere.
 */

#include "MyCollisions.h"
#define alpha 3.0

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
#undef __FUNCT__
#define __FUNCT__ "Vin"

PetscReal Vin(PetscInt i, PetscInt n, PetscReal nn, PetscReal Ti, PetscReal Tn)
{
	PetscReal  Tr,nu;
	
	// Formulas from Schunk and Nagy, 2000, pp. 97--99, Tables 4.4, 4.5, 4.6
	Tr = (Ti+Tn)/2.0;
	if      (i==O2p  && n==CO2) {nu = 5.63e-10*nn*1e-6;}
	else if (i==CO2p && n==CO2) {nu = (Tr>850)*2.85e-11*(nn*1e-6)*sqrt(Tr)*(1-0.083*log10(Tr))*(1-0.083*log10(Tr));}
	else if (i==Op   && n==CO2) {nu = 8.95e-10*nn*1e-6;}
	else if (i==O2p  && n==O)   {nu = 2.31e-10*nn*1e-6;}
	else if (i==CO2p && n==O)   {nu = 1.76e-10*nn*1e-6;}
	else if (i==Op   && n==O)   {nu = (Tr>235)*3.67e-11*(nn*1e-6)*sqrt(Tr)*(1-0.064*log10(Tr))*(1-0.064*log10(Tr));}
	else    {printf("Warning: unknown collision\n\tAssuming null collision frequency.\n"); nu=0.0;}
	
	if(!(nu>=0)) nu = 0.0;
	
	return nu;
}

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
#undef __FUNCT__
#define __FUNCT__ "vin"

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
#undef __FUNCT__
#define __FUNCT__ "Ven"
PetscReal Ven(PetscInt n, PetscReal nn, PetscReal Te)
{
	PetscReal  nu;
	
	// Formulas from Schunk and Nagy, 2000, pp. 97--99, Tables 4.4, 4.5, 4.6
	if      (n==CO2) {nu = 3.68e-8*(nn*1e-6)*(1+4.1e-11*pow(fabs(4500 - Te),2.93));}
	else if (n==O)   {nu = 8.9e-11*(nn*1e-6)*(1+5.7e-4*Te)*sqrt(Te);}
	else    {printf("Warning: unknown collision\n\tAssuming null collision frequency.\n"); nu=0.0;}
	
	if(!(nu>=0)) nu = 0.0;
	
	return nu;
}

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
#undef __FUNCT__
#define __FUNCT__ "ven"

PetscReal ven(PetscReal profiles[][Nz_REF],PetscInt n, PetscReal h)
{
	PetscReal nn=0.0,Te;
	
	Te = Interpolate(profiles,8,h,lin_flat);
	
	if (n==CO2) nn = Interpolate(profiles,5,h,lin_exp);
	if (n==O)   nn = Interpolate(profiles,4,h,lin_exp);
	
	return Ven(n,nn,Te);
}

/*
 * Collisions O2+ CO2 
 */
#undef __FUNCT__
#define __FUNCT__ "v11"
PetscReal v11(PetscReal nCO2){
	PetscReal v;
	nCO2 *= 1e-6; // _cm-3
	v     = 5.63*nCO2; // x1e10_s^-1
	return v*1e-10; // 1/_s
}

/*
 * Collisions O2+ O 
 */
#undef __FUNCT__
#define __FUNCT__ "v12"
PetscReal v12(PetscReal nO  ){
	PetscReal v;
	nO *= 1e-6; // _cm-3
	v   = 2.31*nO; // x1e10_s^-1
	return v*1e-10; // 1/_s
}

/*
 * Collisions O2+ O2+ 
 */
#undef __FUNCT__
#define __FUNCT__ "v13"
PetscReal v13(PetscReal nO2p , PetscReal TO2p ){
	PetscReal v;
	nO2p *= 1e-6; // _cm^-3
	v     = 0.16*nO2p/PetscPowScalar(TO2p,3./2.); // _s^-1
	return v; // _s^-1
}

/*
 * Collisions O2+ CO2+
 */
#undef __FUNCT__
#define __FUNCT__ "v14"
PetscReal v14(PetscReal nCO2p, PetscReal TCO2p){
	PetscReal v;
	nCO2p *= 1e-6; // _cm^-3
	v      = 0.17*nCO2p/PetscPowScalar(TCO2p,3./2.); // _s^-1
	return v; // _s^-1
}

/*
 * Collisions O2+ O+
 */
#undef __FUNCT__
#define __FUNCT__ "v15"
PetscReal v15(PetscReal nOp  , PetscReal TOp  ){
	PetscReal v;
	nOp *= 1e-6; // _cm^-3
	v    = 0.13*nOp/PetscPowScalar(TOp,3./2.); // _s^-1
	return v; // _s^-1
}

/*
 * Collisions O2+ e
 */
#undef __FUNCT__
#define __FUNCT__ "v16"
PetscReal v16(){
	return 0;
}

/*
 * Collisions CO2+ CO2
 */
#undef __FUNCT__
#define __FUNCT__ "v21"
PetscReal v21(PetscReal nCO2 , PetscReal TCO2  , PetscReal TCO2p){
	PetscReal v,Tr;
	nCO2 *= 1e-6; // _cm^-3
	Tr = (TCO2+TCO2p)/2.0; // _K
	v  = (Tr>850)*2.85e-11*nCO2*PetscSqrtScalar(Tr)*PetscPowScalar(1-0.083*log10(Tr),2.0); // _s^-1
	return v; // _s^-1
}

/*
 * Collisions CO2+ O 
 */
#undef __FUNCT__
#define __FUNCT__ "v22"
PetscReal v22(PetscReal nO  ){
	PetscReal v;
	nO *= 1e-6; // _cm-3
	v   = 1.76*nO; // x1e10_s^-1
	return v*1e-10; // 1/_s
}

/*
 * Collisions CO2+ O2+ 
 */
#undef __FUNCT__
#define __FUNCT__ "v23"
PetscReal v23(PetscReal nO2p , PetscReal TO2p ){
	PetscReal v;
	nO2p *= 1e-6; // _cm^-3
	v     = 0.12*nO2p/PetscPowScalar(TO2p,3./2.); // _s^-1
	return v; // _s^-1
}

/*
 * Collisions CO2+ CO2+
 */
#undef __FUNCT__
#define __FUNCT__ "v24"
PetscReal v24(PetscReal nCO2p, PetscReal TCO2p){
	PetscReal v;
	nCO2p *= 1e-6; // _cm^-3
	v      = 0.14*nCO2p/PetscPowScalar(TCO2p,3./2.); // _s^-1
	return v; // _s^-1
}

/*
 * Collisions CO2+ O+
 */
#undef __FUNCT__
#define __FUNCT__ "v25"
PetscReal v25(PetscReal nOp  , PetscReal TOp  ){
	PetscReal v;
	nOp *= 1e-6; // _cm^-3
	v    = 0.10*nOp/PetscPowScalar(TOp,3./2.); // _s^-1
	return v; // _s^-1
}

/*
 * Collisions CO2+ e
 */
#undef __FUNCT__
#define __FUNCT__ "v26"
PetscReal v26(){
	return 0;
}

/*
 * Collisions O+ CO2
 */
#undef __FUNCT__
#define __FUNCT__ "v31"
PetscReal v31(PetscReal nCO2){
	PetscReal v;
	nCO2 *= 1e-6; // _cm^-3
	v     = 8.95*nCO2; // x1e10 _s^-1
	return v*1e-10; // _s^-1
}

/*
 * Collisions O+ O 
 */
#undef __FUNCT__
#define __FUNCT__ "v32"
PetscReal v32(PetscReal nO   , PetscReal TO    , PetscReal TOp  ){
	PetscReal v,Tr;
	nO *= 1e-6; // _cm-3
	Tr  = (TO+TOp)/2.0; // _K
	v   = (Tr>235)*3.67e-11*nO*PetscSqrtScalar(Tr)*PetscPowScalar(1-0.064*log10(Tr),2.0); // _s^-1
	return v; // 1/_s
}

/*
 * Collisions O+ O2+ 
 */
#undef __FUNCT__
#define __FUNCT__ "v33"
PetscReal v33(PetscReal nO2p , PetscReal TO2p ){
	PetscReal v;
	nO2p *= 1e-6; // _cm^-3
	v     = 0.26*nO2p/PetscPowScalar(TO2p,3./2.); // _s^-1
	return v; // _s^-1
}

/*
 * Collisions O+ CO2+
 */
#undef __FUNCT__
#define __FUNCT__ "v34"
PetscReal v34(PetscReal nCO2p, PetscReal TCO2p){
	PetscReal v;
	nCO2p *= 1e-6; // _cm^-3
	v      = 0.27*nCO2p/PetscPowScalar(TCO2p,3./2.); // _s^-1
	return v; // _s^-1
}

/*
 * Collisions O+ O+
 */
#undef __FUNCT__
#define __FUNCT__ "v35"
PetscReal v35(PetscReal nOp  , PetscReal TOp  ){
	PetscReal v;
	nOp *= 1e-6; // _cm^-3
	v    = 0.22*nOp/PetscPowScalar(TOp,3./2.); // _s^-1
	return v; // _s^-1
}

/*
 * Collisions O+ e
 */
#undef __FUNCT__
#define __FUNCT__ "v36"
PetscReal v36(){
	return 0;
}

/*
 * Collisions e CO2
 */
#undef __FUNCT__
#define __FUNCT__ "v41"
PetscReal v41(PetscReal nCO2, PetscReal Te){
	PetscReal v;
	nCO2 *= 1e-6; // _cm^-3
	v     = 3.68e-8*nCO2*(1+4.1e-11*PetscPowScalar(PetscAbsScalar(4500.-Te),2.93)); // _s^-1
	return v; // _s^-1
}

/*
 * Collisions e O 
 */
#undef __FUNCT__
#define __FUNCT__ "v42"
PetscReal v42(PetscReal nO   , PetscReal Te  ){
	PetscReal v;
	nO *= 1e-6; // _cm-3
	v   = 8.9e-11*nO*(1+5.7e-4*Te)*PetscSqrtScalar(Te); // _s^-1
	return v; // 1/_s
}

/*
 * Collisions e O2+ 
 */
#undef __FUNCT__
#define __FUNCT__ "v43"
PetscReal v43(PetscReal nO2p , PetscReal Te   ){
	PetscReal v;
	nO2p *= 1e-6; // _cm^-3
	v     = 54.5*nO2p/PetscPowScalar(Te,3./2.); // _s^-1
	return v; // _s^-1
}

/*
 * Collisions e CO2+
 */
#undef __FUNCT__
#define __FUNCT__ "v44"
PetscReal v44(PetscReal nCO2p, PetscReal Te   ){
	PetscReal v;
	nCO2p *= 1e-6; // _cm^-3
	v      = 54.5*nCO2p/PetscPowScalar(Te,3./2.); // _s^-1
	return v; // _s^-1
}

/*
 * Collisions e O+
 */
#undef __FUNCT__
#define __FUNCT__ "v45"
PetscReal v45(PetscReal nOp  , PetscReal Te   ){
	PetscReal v;
	nOp *= 1e-6; // _cm^-3
	v    = 54.5*nOp/PetscPowScalar(Te,3./2.); // _s^-1
	return v; // _s^-1
}

/*
 * Collisions e e
 */
#undef __FUNCT__
#define __FUNCT__ "v46"
PetscReal v46(PetscReal ne   , PetscReal Te   ){
	PetscReal v;
	ne *= 1e-6; // _cm^-3
	v   = 54.5/PetscSqrtScalar(2.0)*ne/PetscPowScalar(Te,3./2.); // _s^-1
	return v; // _s^-1
	return 0;
}
