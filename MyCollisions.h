/*
 * MyCollisions.h
 * Created by Jeremy A Riousset on 10/30/13
 *
 * Define collision frequencies between ion, electron, and neutral particles in the Martian ionosphere.
 *
 */

#ifndef MYCOLLISIONS_H
#define MYCOLLISIONS_H

#include "MyProfiles.h"

extern PetscReal vin(PetscReal profiles[][Nz_REF],PetscInt i, PetscInt n, PetscReal h);
extern PetscReal Vin(PetscInt i, PetscInt n, PetscReal nn, PetscReal Ti, PetscReal Tn);
extern PetscReal ven(PetscReal profiles[][Nz_REF],PetscInt n, PetscReal h);
extern PetscReal Ven(PetscInt n, PetscReal nn, PetscReal Te);

extern PetscReal v11(PetscReal nCO2);
extern PetscReal v12(PetscReal nO  );
extern PetscReal v13(PetscReal nO2p , PetscReal TO2p );
extern PetscReal v14(PetscReal nCO2p, PetscReal TCO2p);
extern PetscReal v15(PetscReal nOp  , PetscReal TOp  );
extern PetscReal v16();

extern PetscReal v21(PetscReal nCO2 , PetscReal TCO2  , PetscReal TCO2p);
extern PetscReal v22(PetscReal nO  );
extern PetscReal v23(PetscReal nO2p , PetscReal TO2p );
extern PetscReal v24(PetscReal nCO2p, PetscReal TCO2p);
extern PetscReal v25(PetscReal nOp  , PetscReal TOp  );
extern PetscReal v26();

extern PetscReal v31(PetscReal nCO2);
extern PetscReal v32(PetscReal nO   , PetscReal TO    , PetscReal TOp  );
extern PetscReal v33(PetscReal nO2p , PetscReal TO2p );
extern PetscReal v34(PetscReal nCO2p, PetscReal TCO2p);
extern PetscReal v35(PetscReal nOp  , PetscReal TOp  );
extern PetscReal v36();

extern PetscReal v41(PetscReal nCO2 , PetscReal Te   );
extern PetscReal v42(PetscReal nO   , PetscReal Te   );
extern PetscReal v43(PetscReal nO2p , PetscReal Te   );
extern PetscReal v44(PetscReal nCO2p, PetscReal Te   );
extern PetscReal v45(PetscReal nOp  , PetscReal Te   );
extern PetscReal v46(PetscReal ne   , PetscReal Te   );

#endif
