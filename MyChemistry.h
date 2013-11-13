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

extern PetscReal v1();
extern PetscReal k2(PetscReal Te);
extern PetscReal v3();
extern PetscReal k4(PetscReal Te);
extern PetscReal k5(PetscReal Te);
extern PetscReal k6(PetscReal Te);
extern PetscReal k7(PetscReal Te);
extern PetscReal k8();
extern PetscReal k9();
extern PetscReal k10();

#endif
