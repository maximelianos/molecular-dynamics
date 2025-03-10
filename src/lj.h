/*
* Copyright 2024 Maxim Velikanov
* MIT license
*
* Compute the forces using Lennard-Jones potential.
*
*/

#ifndef LJ_H
#define LJ_H

#include <iostream>

#include "verlet.h"
#include "neighbors.h"

double lj_neighbor_list(Atoms &atoms, NeighborList &neighbor_list, double cutoff, double epsilon = 1.0, double sigma = 1.0);

#endif //LJ_H
