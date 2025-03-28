//
// Created by Maksim Velikanov on 10/07/24.
//

#include "verlet.h"
#include <iostream>

void verlet_step1(Atoms &atoms, double timestep, double mass) {
    atoms.velocities += 0.5 * atoms.forces * timestep / mass;
    atoms.positions += atoms.velocities * timestep;
}

void verlet_step2(Atoms &atoms, double timestep, double mass) {
    // recompute forces before this!
    atoms.velocities += 0.5 * atoms.forces * timestep / mass;
}

