//
// Created by iter on 11.09.24.
//

#ifndef THERMOSTAT_H
#define THERMOSTAT_H

#include "atoms.h"

void berendsen_thermostat(Atoms &atoms, double temperature, double timestep,
                          double relaxation_time);


#endif //THERMOSTAT_H
