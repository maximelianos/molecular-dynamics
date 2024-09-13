//
// Copyright 2024 Maksim Velikanov
//

#ifndef THERMOSTAT_H
#define THERMOSTAT_H

#include "atoms.h"

// adjust atoms.velocities to prevent explosion
void berendsen_thermostat(Atoms &atoms, double temperature, double timestep,
                          double relaxation_time, double mass=1.0);


#endif //THERMOSTAT_H
