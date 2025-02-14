//
// Created by iter on 05.09.24.
//

#ifndef LJ_DIRECT_SUMMATION_H
#define LJ_DIRECT_SUMMATION_H

#include <iostream>

#include "verlet.h"
#include "lj_direct_summation.h"

double kinetic_energy(Atoms &atoms, double mass=1.0);
double get_temperature(Atoms &atoms);
double lj_direct_summation(Atoms &atoms, double epsilon = 1.0, double sigma = 1.0);

#endif //LJ_DIRECT_SUMMATION_H
