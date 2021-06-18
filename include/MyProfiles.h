/*
 * MyProfiles.h
 * Created by Jeremy A Riousset & Carol S Paty on 03/16/11
 *
 * Define a set of very useful profiles between 50 and 450 km.
 *
 * The  altitude grid comes from the transport model, in kilometers,
 * starting from 400 down to 50.
 * The profiles include:
 * -Number density of neutral atomic oxygen in cm^-3 from
 *      Bougher's model (equinox, solar moderate conditions, 2.5 degrees N
 *      latitude, 03 local time).
 * - Density of CO2: same as number density but for CO2.
 * - Total density assuming O and CO2 only (i.e.,density_o + density_co2).
 * - Ionization_rate: this is actually the average ionization rate for
 *      the above atmosphere for 100 different precipitating electron energy
 *      spectra. There is probably nowhere on the planet with this exact
 *      ionization rate profile, it is just the average of Matt's work from the SWIM conference.
 * - Electron density in cm^-3 from the production rate
 *      above (assuming all ions are O2+ and using the O2+ temperature dependent
 *      dissociative recombination rate).
 * - Neutral temperature in Kelvin from Bougher
 * - Ion temperature in Kelvin. Matt made it up to roughly match the ion
 *      temperature profile from Fox, 1993 (more or less). Complete fabrication
 *      based loosely on the data.
 * - Electron temperature in Kelvin. Matt just used the
 *      equations given in Fox, 1993.
 * - x, y, and z components of the magnetic
 *      field in nT. Similar to ionization_rate.txt, these are the average
 *      profiles from 100 different (but close) locations. The individual
 *      profiles were taken from the Cain model (with n = 60). Probably nowhere
 *      on Mars has this exact profile, but it's representative of a large area.
 * - Magnitude of the average magnetic field from above
 *      (not the average of the magnitudes). It is a little bit weaker than the
 *      1000 nT field we estimated (only 370 nT), but these are the fields (and
 *      atmosphere) that Matt used to calculate currents for his SWIM stuff, so we
 *      should still get a nice dynamo region.
 */

#ifndef MYPROFILES_H
#define MYPROFILES_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "petscsys.h"
#include <assert.h>
#include "MyCtx.h"

/*
 * To extend the model to N-ion species, we shall proceed as follows:
 * 1. remove "e" from the charged species enum list
 * 2. correct all the errors
 * 3. rewrite the model to include the equation for electron distincly from the ion equations
 * 4. add the new ion(s) to the chargedspecies enum list
 * 5. modify the input files accordingly
 *    - Add a line in main.in to enter the number of ion species
 *    - Add the adequate number of masses in main.in
 *    - Automatize the creation of the DMDA da to use the appropriate number of equations and variables
 * 6. Modify the size of the containers ni, pi, vi... in MyCtx.h
 */
enum interpolations {lin_lin, lin_exp, lin_flat};

extern PetscErrorCode ReadTable(PetscReal table[][Nz_REF], PetscInt P, const char *fName);
extern PetscReal Interpolate(PetscReal table[][Nz_REF], PetscInt i, PetscReal h, PetscInt ItpType);  
extern PetscReal Interpolate1(PetscReal x, PetscReal xi[], PetscReal yi[], PetscInt Ni, PetscInt ItpType);  
extern PetscReal Interpolate2(PetscReal x, PetscReal xi[], PetscReal yi[], PetscInt Ni, PetscInt ItpType);  

#endif

