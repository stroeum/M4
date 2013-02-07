/*
 * MyPDEs.h
 * Created by Jeremy A Riousset on 02/04/11
 * Define PDEs system
 */

#ifndef MYPDES_H
#define MYPDES_H

#include "MyCtx.h"

extern PetscErrorCode FormBCu(Vec,void*);
extern PetscErrorCode CheckValues(Vec,Vec,void*);
extern PetscErrorCode FormIntermediateFunction(Vec,Vec,void*);
extern PetscErrorCode FormFunction(TS,PetscReal,Vec,Vec,void*);  
extern PetscErrorCode FormInitialSolution(Vec,void*);

#endif

