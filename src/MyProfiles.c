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
PetscErrorCode ReadTable(PetscReal table[][Nz_REF], PetscInt P, const char *fName)
{
	char buffer[MAX_LINE_LENGTH];
	FILE *pFile;
	int m=0, p=0;

	PetscFunctionBegin;
	
	pFile=fopen(fName,"r");
	fgets(buffer, MAX_LINE_LENGTH, pFile);
	
	for (m=0; m<Nz_REF; m++) {
		for (p=0; p<P; p++) {
			fscanf(pFile,"%lf", &table[p][m]);
		}
	}
	fclose(pFile);
	
 	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "Interpolate1" // xi is in ascending order
PetscReal Interpolate1(PetscReal x, PetscReal xi[], PetscReal yi[], PetscInt Ni, PetscInt ItpType) {  
	PetscInt id[3];
	PetscReal y=0.0, L=0.0, dx=0.0, dl[2]={0.0,0.0}, Dx[3]={0.0,0.0,0.0}, Dy=0.0;
	PetscReal X[3]={0.0,0.0,0.0},Y[3]={0.0,0.0,0.0};
	
	if (x>=xi[0] && x<=xi[Ni-1]) {
		id[0]=0;
		id[2]=Ni-1;
		while(id[2]!=id[0]+1) {
			id[1] = (PetscInt)((id[0]+id[2])/2);
			dx = x-xi[id[1]];
			if(dx <0) id[2]=id[1];
			if(dx>=0) id[0]=id[1];
		}
		X[0] = xi[id[0]];
		X[1] = xi[id[2]];
		Y[0] = yi[id[0]];
		Y[1] = yi[id[2]];
		
		y = (x-X[0])/(X[1]-X[0])*Y[1] + (x-X[1])/(X[0]-X[1])*Y[0];
	} else if (x>xi[Ni-1] || x<xi[0]) {
		if (x>xi[Ni-1]) {
			id[0] = Ni-1;
			id[1] = Ni-2;
			id[2] = Ni-3;
		} else {
			id[0] = 0;
			id[1] = 1;
			id[2] = 2;
		}
		X[0]  = xi[id[0]]; X[1]  = xi[id[1]]; X[2]  = xi[id[2]];
		Y[0]  = yi[id[0]]; Y[1]  = yi[id[1]]; Y[2]  = yi[id[2]];
		dl[0] = X[1] - X[0];
		dl[1] = X[2] - X[0];
		Dx[0] = - 1.0/dl[0] - 1.0/dl[1]; 
		Dx[1] =   1.0/dl[0] - 1.0/(dl[0]-dl[1]); 
		Dx[2] =   1.0/dl[1] + 1.0/(dl[0]-dl[1]);
		
		Dy    = Dx[0]*Y[0] + Dx[1]*Y[1] + Dx[2]*Y[2];
		
		// linear interpolation above max xi
		if (ItpType==lin_lin) {
			y = Dy*(x-X[0]) + Y[0];
		}
		// exponential decrease above max xi
		if (ItpType==lin_exp) {
			L = Y[0]/Dy;
			if (L<0) {
				if (fabs((x-X[0])/L)<alpha)
					y = Y[0]*exp( (x-X[0])/L );
				else
					y = Y[0]*exp(-alpha);
			} else {
				PetscPrintf(PETSC_COMM_WORLD,"Wrong Interpolation\n");
				exit(1);
			}
		}
		// flat interpolation above max xi
		if (ItpType==lin_flat) {
			y = Y[0];
		}
	} else {
		printf("h is outside of interpolation range\n");
		exit(0);
	}
	return y;
}

#undef __FUNCT__
#define __FUNCT__ "Interpolate2" // xi is in descending order
PetscReal Interpolate2(PetscReal x, PetscReal xi[], PetscReal yi[], PetscInt Ni, PetscInt ItpType) {  
	PetscInt id[3];
	PetscReal y=0.0, L=0.0, dx=0.0, dl[2]={0.0,0.0}, Dx[3]={0.0,0.0,0.0}, Dy=0.0;
	PetscReal X[3]={0.0,0.0,0.0},Y[3]={0.0,0.0,0.0};
	
	if (x>=xi[Ni-1] && x<=xi[0]) {
		id[2]=0;
		id[0]=Ni-1;
		while(id[2]!=id[0]-1) {
			id[1] = (PetscInt)((id[0]+id[2])/2);
			dx = x-xi[id[1]];
			if(dx <0) id[2]=id[1];
			if(dx>=0) id[0]=id[1];
		}
		X[0] = xi[id[0]];
		X[1] = xi[id[2]];
		Y[0] = yi[id[0]];
		Y[1] = yi[id[2]];
		
		y = (x-X[0])/(X[1]-X[0])*Y[1] + (x-X[1])/(X[0]-X[1])*Y[0];
	} else if (x>xi[0]) {
		id[0] = 0;
		id[1] = 1;
		id[2] = 2;
		X[0]  = xi[id[0]]; X[1]  = xi[id[1]]; X[2]  = xi[id[2]];
		Y[0]  = yi[id[0]]; Y[1]  = yi[id[1]]; Y[2]  = yi[id[2]];
		dl[0] = X[1] - X[0];
		dl[1] = X[2] - X[0];
		Dx[0] = - 1.0/dl[0] - 1.0/dl[1]; 
		Dx[1] =   1.0/dl[0] - 1.0/(dl[0]-dl[1]); 
		Dx[2] =   1.0/dl[1] + 1.0/(dl[0]-dl[1]);
		
		Dy    = Dx[0]*Y[0] + Dx[1]*Y[1] + Dx[2]*Y[2];
		
		// linear interpolation above max reference altitude
		if (ItpType==lin_lin) {
			y = Dy*(x-X[0]) + Y[0];
		}
		// exponential decrease above max reference altitude
		if (ItpType==lin_exp) {
			L = Y[0]/Dy;
			if (L<0) {
				if (fabs((x-X[0])/L)<alpha)
					y = Y[0]*exp( (x-X[0])/L );
				else
					y = Y[0]*exp(-alpha);
			} else {
				PetscPrintf(PETSC_COMM_WORLD,"Wrong Interpolation\n");
				exit(1);
			}
		}
		// flat interpolation above max reference altitude
		if (ItpType==lin_flat) {
			y = Y[0];
		}
	} else if (x<xi[Ni-1]) {
		id[0] = Ni-1;
		id[1] = Ni-2;
		id[2] = Ni-3;
		X[0]  = xi[id[0]]; X[1]  = xi[id[1]]; X[2]  = xi[id[2]];
		Y[0]  = yi[id[0]]; Y[1]  = yi[id[1]]; Y[2]  = yi[id[2]];
		dl[0] = X[1] - X[0];
		dl[1] = X[2] - X[0];
		Dx[0] = - 1.0/dl[0] - 1.0/dl[1]; 
		Dx[1] =   1.0/dl[0] - 1.0/(dl[0]-dl[1]); 
		Dx[2] =   1.0/dl[1] + 1.0/(dl[0]-dl[1]);
		
		Dy    = Dx[0]*Y[0] + Dx[1]*Y[1] + Dx[2]*Y[2];
		
		// linear interpolation below min reference altitude
		if (ItpType==lin_lin) {
			y = Dy*(x-X[0]) + Y[0];
		}
		// exponential decrease below min reference altitude
		if (ItpType==lin_exp) {
			L = Y[0]/Dy;
			if (fabs((x-X[0])/L)<alpha)
				y = Y[0]*exp( (x-X[0])/L );
			else
				y = Y[0]*exp(-alpha);
		}
		// flat interpolation below min reference altitude
		if (ItpType==lin_flat) {
			y = Y[0];
		}
	} else {
		printf("h is outside of interpolation range\n");
		exit(0);
	}
	return y;
}

#undef __FUNCT__
#define __FUNCT__ "Interpolate"
PetscReal Interpolate(PetscReal profiles[][Nz_REF], PetscInt i, PetscReal h, PetscInt ItpType)
{
	return Interpolate2(h,profiles[0],profiles[i+1],Nz_REF,ItpType);
}

