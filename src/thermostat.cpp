//
// Copyright 2024 Maksim Velikanov
//

#include "thermostat.h"

#include <lj_direct_summation.h>

void berendsen_thermostat(Atoms &atoms, double temperature, double target_temperature, double timestep,
                          double relaxation_time, double mass) {
    double lambda = std::sqrt(1 + (target_temperature / temperature - 1) * timestep / relaxation_time);
    atoms.velocities = atoms.velocities * lambda;
}