//
// Created by iter on 05.09.24.
//

#ifndef LJ_H
#define LJ_H

#include <iostream>

#include "verlet.h"
#include "neighbors.h"

double lj_neighbor_list(Atoms &atoms, NeighborList &neighbor_list, double epsilon = 1.0, double sigma = 1.0);

#endif //LJ_H
