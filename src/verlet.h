/*
* Copyright 2024 Maxim Velikanov
* MIT license
*
* Verlet integrator recomputes the positions of atoms.
*
*/

#ifndef MD_MASTER_VERLET_H
#define MD_MASTER_VERLET_H

#include "atoms.h"

void verlet_step1(Atoms &atoms, double timestep, double mass);
void verlet_step2(Atoms &atoms, double timestep, double mass);

void run_simulation(Atoms &atoms, int nb_steps, double timestep, double mass);


#endif // MD_MASTER_VERLET_H
