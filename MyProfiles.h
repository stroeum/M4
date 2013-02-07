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

enum chargedspecies {O2p, CO2p, Op, e};
enum neutrals {CO2, O};

extern PetscReal Profile(int,double);  
extern PetscReal Partition(int,double);  
extern PetscReal vin(PetscInt i, PetscInt n, PetscReal h);
extern PetscReal Vin(PetscInt i, PetscInt n, PetscReal nn, PetscReal Ti, PetscReal Tn);
extern PetscReal ven(PetscInt n, PetscReal h);
extern PetscReal Ven(PetscInt n, PetscReal nn, PetscReal Te);

#endif

