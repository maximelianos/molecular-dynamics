//
// Created by Maksim Velikanov on 10/07/24.
//

#ifndef MD_MASTER_VERLET_H
#define MD_MASTER_VERLET_H

#include "atoms.h"

void verlet_step1(Atoms &atoms, double timestep, double mass);
void verlet_step2(Atoms &atoms, double timestep, double mass);

void run_simulation(Atoms &atoms, int nb_steps, double timestep, double mass);


#endif // MD_MASTER_VERLET_H
