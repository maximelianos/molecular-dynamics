//
// Created by iter on 05.09.24.
//

#include "verlet.h"
#include "lj_direct_summation.h"
#include <iostream>
#include <gtest/gtest.h>

TEST(LennardJonesTest, Norm) {
    Point_t a;
    a << 1, 0, 0;
    Point_t b;
    b << 0, 1, 0;
    ASSERT_FLOAT_EQ((a - b).norm(), std::sqrt(2));
}

// Demonstrate some basic assertions.
TEST(LennardJonesTest, TwoAtoms) {
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

    lj_direct_summation(atoms);

    std::cout << atoms.forces;

    //ASSERT_FLOAT_EQ(atoms.positions(1, 1), -24);
}
