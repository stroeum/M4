/*
 * MyMonitor.h
 * Created by Jeremy A Riousset on 02/04/11
 * Output solutions to binary files for with PDE solvers
 */

#ifndef MYMONITOR_H
#define MYMONITOR_H

#include "MyCtx.h"
#include "MyPDEs.h"

extern PetscErrorCode MyTSMonitor(TS,PetscInt,PetscReal,Vec,void*);
extern PetscErrorCode MySNESMonitor(SNES,PetscInt,PetscReal,void*);

#endif
