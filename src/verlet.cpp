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
    atoms.velocities += 0.5 * atoms.forces * timestep / mass;
}

void run_simulation(Atoms &atoms, int nb_steps) {
    // s = v_0 * t + (a * t^2) / 2
    for (int i = 0; i < nb_steps; i++) {
        verlet_step1(atoms, 1, 1);
        verlet_step2(atoms, 1, 1);
    }
}


