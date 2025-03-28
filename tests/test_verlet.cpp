// Include a library file to make sure proper includes are set
#include "hello.h"
#include "verlet.h"
#include <iostream>
#include <gtest/gtest.h>

void run_simulation(Atoms &atoms, int nb_steps, double timestep=1.0, double mass=1.0) {
    for (int i = 0; i < nb_steps; i++) {
        verlet_step1(atoms, timestep, mass);
        verlet_step2(atoms, timestep, mass);
    }
}

// Demonstrate some basic assertions.
TEST(VerletTest, ConstantForce) {
    // atom 1
    Positions_t pos(3, 2);
    pos(0, 0) = 1; // x
    pos(1, 0) = 0; // y
    pos(2, 0) = 0; // z
    // atom 2
    pos(0, 1) = 0;
    pos(1, 1) = 1;
    pos(2, 1) = 0;

    Atoms atoms(pos);
    // atom 1
    atoms.forces(0, 0) = 1;
    atoms.forces(1, 0) = 0;
    atoms.forces(2, 0) = 0;
    // atom 2
    atoms.forces(0, 1) = 0;
    atoms.forces(1, 1) = -0.5;
    atoms.forces(2, 1) = 0;

    run_simulation(atoms, 10, 1.0, 1.0);

    // s = v_0 * t + (a * t^2) / 2
    // atom 1
    // x_1 = x_0 + s =
    // 1 + 0 + (1 * 100) / 2 = 51
    ASSERT_FLOAT_EQ(atoms.positions(0, 0), 51);
    // atom 2
    // f = ma, a = f/m
    // y_1 = y_0 + s =
    // 1 + 0 + (-0.5 * 100) / 2 = 1 - 25 = 25
    ASSERT_FLOAT_EQ(atoms.positions(1, 1), -24);
}


TEST(VerletTest, LinearForce) {
    // Linearly changing force F(t) = 1 - t
    // S(t) = t^2 / 2 - t^3 / 6, when v_0 = 0
    // atom 1
    Positions_t pos(3, 2);
    pos(0, 0) = 1; // x
    pos(1, 0) = 0; // y
    pos(2, 0) = 0; // z

    Atoms atoms(pos);
    // atom 1
    atoms.forces(0, 0) = 0;
    atoms.forces(1, 0) = 0;
    atoms.forces(2, 0) = 0;

    // end t = 0.5
    // S(0.5) = 5 / (8*6)

    int nb_steps = 50000;
    double timestep = 0.00001;
    double mass = 1.0;
    double t = 0;
    for (int i = 0; i < nb_steps; i++) {
        atoms.forces(0, 0) = 1 - t;
        verlet_step1(atoms, timestep, mass);
        t = timestep * (i + 0.5);

        atoms.forces(0, 0) = 1 - t;
        verlet_step2(atoms, timestep, mass);
        t = timestep * (i + 1);
    }

    // atom 1
    ASSERT_FLOAT_EQ(atoms.positions(0, 0), 1 + 5.0/(8.0*6.0));
}
