// Include a library file to make sure proper includes are set
#include "hello.h"
#include "thermostat.h"
#include "verlet.h"
#include <gtest/gtest.h>
#include <iostream>
#include <lj_direct_summation.h>

// Demonstrate some basic assertions.
TEST(ThermostatTest, OneAtom) {
    // atom 1
    Positions_t pos(3, 2);
    pos(0, 0) = 1; // x
    pos(1, 0) = 0; // y
    pos(2, 0) = 0; // z

    Atoms atoms(pos);
    // atom 1
    atoms.velocities(0, 0) = 1;
    atoms.velocities(1, 0) = 0;
    atoms.velocities(2, 0) = 0;

    double sigma = 1;
    double epsilon = 1;
    double m = 1;

    double step_t = 1 * std::sqrt(m * sigma * sigma / epsilon);

    // velocity rescale - when relaxation time = step time
    berendsen_thermostat(atoms, 0.5, step_t, step_t);
    EXPECT_FLOAT_EQ(get_temperature(atoms), 0.5);

    int n_steps = 50;
    // achieve temperature 1000.0
    for (int i = 0; i < n_steps; i++) {
        // apply forces
        verlet_step1(atoms, step_t, m);
        verlet_step2(atoms, step_t, m);
        berendsen_thermostat(atoms, 1000.0, step_t, step_t * 10);
    }
    EXPECT_NEAR(get_temperature(atoms) / 1000.0, 1.0, 1e-2);
}
