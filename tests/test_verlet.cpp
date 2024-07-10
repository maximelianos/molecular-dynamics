// Include a library file to make sure proper includes are set
#include "hello.h"
#include "verlet.h"
#include <iostream>
#include <gtest/gtest.h>

// Demonstrate some basic assertions.
TEST(VerletTest, ConstantForce) {
    Positions_t pos(3, 1);
    pos(0, 0) = 1;
    pos(1, 0) = 0;
    pos(2, 0) = 0;
    Atoms atoms(pos);
    atoms.forces(0, 0) = 1;
    atoms.forces(1, 0) = 0;
    atoms.forces(2, 0) = 0;

    int nb_steps = 10;
    for (int i = 0; i < nb_steps; i++) {
        verlet_step1(atoms, 1, 1);
        verlet_step2(atoms, 1, 1);
    }

    ASSERT_FLOAT_EQ(atoms.positions(0, 0), 51);
}
