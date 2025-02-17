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

struct Average {
  double sum = 0;
  int n = 0;

  void add(double value) {
    sum += value;
    n += 1;
  }

  double result() {
    double res = sum / n;
    sum = 0;
    n = 0;
    return res;
  }
};

double kinetic_energy(Atoms &atoms, double mass=1.0);
double get_temperature_lj(double e_kin, double nb_atoms);
double get_temperature(double e_kin, double nb_atoms);
double lj_direct_summation(Atoms &atoms, double epsilon = 1.0, double sigma = 1.0);

#endif //LJ_DIRECT_SUMMATION_H
