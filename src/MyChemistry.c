/*
 * MyChemistry.c
 * Created by Jeremy A Riousset & Carol S Paty on 10/31/13
 *
 * Define chemical reactions in the Martian ionosphere
 *
 */

#include "MyChemistry.h"

/*
 * Photoionization of CO2
 * CO2 + hv -> CO2+ + e
 */
#undef __FUNCT__
#define __FUNCT__ "v1"
PetscReal v1(void) {
	return 5.0e-7; //_s-1
}

/*
 * Electron impact ionization of CO2
 * CO2 + e -> CO2+ + e + e
 */
#undef __FUNCT__
#define __FUNCT__ "v2"
PetscReal v2(PetscReal N , PetscReal Te) {
		// Ei in reverse order //
	/*
	 PetscInt  Ni=58;
	 PetscReal me=9.1093e-31; //_kg
	 PetscReal E=kB*Te/qe;
	 PetscReal Xsec;
	 PetscReal K,v;
	 PetscReal Ei[] = {1000,950,900,850,800,750,700,650,600,550,500,450,400,350,300,275,250,225,200,180,160,140,120,110,100,95,90,85,80,75,70,65,60,55,50,45,40,35,30,25,24,23.5,23,22.5,22,21.5,21,19.5,19,18.5,18,17.5,17,16.5,16,15.5,15,14.5}; //_eV, Electron energy
	 PetscReal XsecCO2i[] = {0.876,0.909,0.941,0.964,1.01,1.06,1.1,1.16,1.21,1.27,1.35,1.43,1.54,1.62,1.75,1.83,1.87,1.95,2.01,2.08,2.12,2.19,2.23,2.23,2.25,2.23,2.22,2.2,2.19,2.15,2.13,2.1,2.06,2,1.94,1.84,1.7,1.53,1.34,0.969,0.88,0.828,0.777,0.727,0.676,0.623,0.577,0.452,0.428,0.373,0.333,0.293,0.255,0.215,0.174,0.135,0.097,0.055}; //1e-16_cm2, CO2+ ionization cross section
	 PetscReal Xsectoti[] = {1.3,1.36,1.41,1.45,1.53,1.61,1.68,1.76,1.85,1.96,2.09,2.23,2.43,2.58,2.82,2.97,3.05,3.21,3.32,3.43,3.52,3.63,3.66,3.66,3.64,3.6,3.56,3.51,3.45,3.36,3.27,3.18,3.06,2.88,2.71,2.5,2.25,1.95,1.58,1.04,0.88,0.828,0.777,0.727,0.676,0.623,0.577,0.452,0.428,0.373,0.333,0.293,0.255,0.215,0.174,0.135,0.097,0.055}; //1e-16_cm2, total ionization cross section
	 Xsec = Interpolate2(E,Ei,XsecCO2i,Ni,lin_lin)*1e-20; //_m2
	 v = PetscSqrtScalar(2*kB*Te/me); //_m/_s
	 K = Xsec*v; //_m3/_s
	 return K;
	 */
	
		// Ei in direct order //
	PetscInt  Ni=58;
	PetscReal me=9.1093e-31; //_kg
    PetscReal E=kB*Te/qe;
//    PetscReal E = 15;
	PetscReal Xsec;
	PetscReal K,v,V;
	PetscReal Ei[] = {14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,21,21.5,22,22.5,23,23.5,24,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,140,160,180,200,225,250,275,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000}; //_eV, Electron energy
	PetscReal XsecCO2i[] = {0.055,0.097,0.135,0.174,0.215,0.255,0.293,0.333,0.373,0.428,0.452,0.577,0.623,0.676,0.727,0.777,0.828,0.88,0.969,1.34,1.53,1.7,1.84,1.94,2,2.06,2.1,2.13,2.15,2.19,2.2,2.22,2.23,2.25,2.23,2.23,2.19,2.12,2.08,2.01,1.95,1.87,1.83,1.75,1.62,1.54,1.43,1.35,1.27,1.21,1.16,1.1,1.06,1.01,0.964,0.941,0.909,0.876}; //1e-16_cm2, CO2+ ionization cross section
//    PetscReal Xsectoti[] = {0.055,0.097,0.135,0.174,0.215,0.255,0.293,0.333,0.373,0.428,0.452,0.577,0.623,0.676,0.727,0.777,0.828,0.88,1.04,1.58,1.95,2.25,2.5,2.71,2.88,3.06,3.18,3.27,3.36,3.45,3.51,3.56,3.6,3.64,3.66,3.66,3.63,3.52,3.43,3.32,3.21,3.05,2.97,2.82,2.58,2.43,2.23,2.09,1.96,1.85,1.76,1.68,1.61,1.53,1.45,1.41,1.36,1.3}; //1e-16_cm2, total ionization cross section
    /*
        Negative values are coming from the Interpolate function. The function is returning negative numbers because E is very low due to Te, Xsec is then interpolated to a negative number
     */
	Xsec = Interpolate1(E,Ei,XsecCO2i,Ni,lin_lin)*1e-20; //_m2,
	if(Xsec<0) Xsec = 0;
	v = PetscSqrtScalar(2*kB*Te/me); //_m/_s
	K = Xsec*v; //_m3/_s
    V = N*K;
//    PetscPrintf(PETSC_COMM_WORLD, "%e %e %e\n", Te,E,Xsec);
    return V; //_m3/_s
}

/*
 * Photoionization of O
 * O + hv -> O+ + e
 */
#undef __FUNCT__
#define __FUNCT__ "v3"
PetscReal v3(void) {
	return 2.0e-7; //_s-1
}

/*
 * Electron impact ionization of O
 * O + e -> O+ + e + e
 */
#undef __FUNCT__
#define __FUNCT__ "v4"
PetscReal v4(PetscReal N , PetscReal Te) {
	PetscReal K,E,v;
	E = kB*Te/qe;
	K = 9.0e-9*PetscPowScalar(E,0.7)*PetscExpScalar(-13.8/E); //_cm^3/_s
    v = N*K;
	return v*1e-6; //_m3/_s
}

/*
 * Dissociative recombination of O2+
 * O2+ + e -> O + O
 */
#undef __FUNCT__
#define __FUNCT__ "v5"
PetscReal v5(PetscReal N , PetscReal Te) {
	PetscReal K,v;
	if (Te<1200) {
		K = 1.95e-7*PetscPowScalar( 300/Te,0.7 ); //_cm^3/_s
	} else {
		K = 7.38e-8*PetscPowScalar(1200/Te,0.56); //_cm^3/_s
	}
    v = N*K;
    return v*1e-6; //_m3/_s
}

/*
 * Dissociative recombination of CO2+
 * CO2+ + e -> CO + O
 */
#undef __FUNCT__
#define __FUNCT__ "v6"
PetscReal v6(PetscReal N , PetscReal Te) {
	PetscReal K, v;
	K = 3.1e-7*PetscPowScalar(300/Te,0.5); //_cm^3/_s
    v = N*K;
    return v*1e-6; //_m3/_s
}

/*
 * Radiative recombination of O+
 * O+ + e -> O + hv
 */
#undef __FUNCT__
#define __FUNCT__ "v7"
PetscReal v7(PetscReal N , PetscReal Te) {
	PetscReal K, v;
	K = 3.71e-12*PetscPowScalar(250/Te,0.7); //_cm^3/_s
    v = N*K;
    return v*1e-6; //_m3/_s
}

/*
 * Positive ion chemistry
 * CO2+ + O -> O2+ + CO
 */
#undef __FUNCT__
#define __FUNCT__ "v8"
PetscReal v8(PetscReal N) {
	PetscReal K, v;
	K = 1.64e-10; //_cm^3/_s
    v = N*K;
    return v*1e-6; //_m3/_s
}

/*
 * Positive ion chemistry
 * CO2+ + O -> O+ + CO2
 */
#undef __FUNCT__
#define __FUNCT__ "v9"
PetscReal v9(PetscReal N) {
	PetscReal K, v;
	K = 9.6e-11; //_cm^3/_s
    v = N*K;
    return v*1e-6; //_m3/_s
}

/*
 * Positive ion chemistry
 * O+ + CO2 -> O2+ + CO
 */
#undef __FUNCT__
#define __FUNCT__ "v10"
PetscReal v10(PetscReal N) {
	PetscReal K, v;
	K = 1.1e-9; //_cm^3/_s
    v = N*K;
    return v*1e-6; //_m3/_s
}
