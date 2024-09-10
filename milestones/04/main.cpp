/*
* Copyright 2024 Maxim Velikanov
*
* ### MIT license
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in
* all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/

#include <iomanip>
#include <iostream>

#include "verlet.h"
#include "lj_direct_summation.h"
#include "xyz.h"



int main() {
    auto [names, positions, velocities]{read_xyz_with_velocities("lj54.xyz")};
    Atoms atoms(positions);
    atoms.velocities = velocities;

    double sigma = 1;
    double epsilon = 1;
    double m = 1;

    double end_t = 100 * std::sqrt(m * sigma * sigma / epsilon);
    double begin_t = 0;
    double step_t = 0.001 * std::sqrt(m * sigma * sigma / epsilon);
    int steps = 0;

    for (; begin_t < end_t; begin_t += step_t) {
        steps += 1;

        // compute forces
        double e_pot = lj_direct_summation(atoms, epsilon, sigma);
        // apply forces
        verlet_step1(atoms, step_t, m);
        verlet_step2(atoms, step_t, m);
        // compute total energy
        double e = e_pot + kinetic_energy(atoms);
        if (steps % 100 == 0)
            std::cout << "e_tot " << e << "\n";
    }
    return 0;
}
