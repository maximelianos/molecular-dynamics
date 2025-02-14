//
// Copyright 2024 Maksim Velikanov
//

#include "thermostat.h"

#include <lj_direct_summation.h>

void berendsen_thermostat(Atoms &atoms, double temperature, double timestep,
                          double relaxation_time, double mass) {
    double t = get_temperature(atoms);
    double lambda = std::sqrt(1 + (temperature / t - 1) * timestep / relaxation_time);
    atoms.velocities = atoms.velocities * lambda;
}