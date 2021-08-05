/*
 * MyPDEs.h
 * Created by Jeremy A Riousset on 02/04/11
 * Define PDEs system
 */

#ifndef MYPDES_H
#define MYPDES_H

#include "MyCtx.h"
#include "MyProfiles.h"
#include "MyChemistry.h"
#include "MyCollisions.h"

extern PetscErrorCode FormBCu(PetscReal****,void*);
extern PetscErrorCode FormBCv(PetscReal****,void*);
extern PetscErrorCode CheckValues(PetscReal****,PetscReal****,void*);
extern PetscErrorCode CalculateFluxes(PetscReal,PetscInt,PetscReal****,PetscReal****,void*);
extern PetscErrorCode FormIntermediateFunction(PetscReal****,Vec,void*);
extern PetscErrorCode FormFunction(TS,PetscReal,Vec,Vec,void*);  
extern PetscErrorCode FormInitialSolution(Vec,void*);
extern PetscErrorCode FormInitialBField(PetscReal*,void*,PetscInt,PetscInt*,PetscReal,PetscReal,PetscReal,PetscReal,PetscReal,PetscReal,PetscReal);

#endif

