/*
 * MyChemistry.h
 * Created by Jeremy A Riousset on 10/30/13
 *
 * Define reaction rates in the Martian atmosphere.
 *
 */

#ifndef MYCHEMISTRY_H
#define MYCHEMISTRY_H

#include "MyCtx.h"
#include "MyProfiles.h"

//extern PetscReal v1(void);
extern PetscReal v1(PetscReal N, PetscReal Z);

extern PetscReal v2(PetscReal N , PetscReal Te);

//extern PetscReal v3(void);
extern PetscReal v3(PetscReal N, PetscReal Z);

extern PetscReal v4(PetscReal N , PetscReal Te);
extern PetscReal v5(PetscReal N , PetscReal Te);
extern PetscReal v6(PetscReal N , PetscReal Te);
extern PetscReal v7(PetscReal N , PetscReal Te);
extern PetscReal v8(PetscReal N);
extern PetscReal v9(PetscReal N);
extern PetscReal v10(PetscReal N);
extern PetscReal k11(void);
extern PetscReal k12(void);
extern PetscReal k13(void);
extern PetscReal k14(void);
extern PetscReal k15(void);

#endif
