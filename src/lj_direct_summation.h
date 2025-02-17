/*
* Copyright 2024 Maxim Velikanov
* MIT license
*
* Compute the forces from Lennard-Jones potential without neighbour lists.
*
*/

#ifndef LJ_DIRECT_SUMMATION_H
#define LJ_DIRECT_SUMMATION_H

#include "verlet.h"
#include "lj_direct_summation.h"

double kinetic_energy(Atoms &atoms, double mass=1.0);
double get_temperature(double e_kin, double nb_atoms);
double lj_direct_summation(Atoms &atoms, double epsilon = 1.0, double sigma = 1.0);

#endif //LJ_DIRECT_SUMMATION_H
