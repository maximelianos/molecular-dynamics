/*
* Copyright 2024 Maxim Velikanov
* MIT license
*
* Adjust the temperature with Berendsen thermostat.
*
*/

#ifndef THERMOSTAT_H
#define THERMOSTAT_H

#include "atoms.h"

// adjust atoms.velocities to prevent explosion
void berendsen_thermostat(Atoms &atoms, double temperature, double target_temperature, double timestep,
                          double relaxation_time, double mass=1.0);


#endif //THERMOSTAT_H
