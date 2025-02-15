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

#include "lj_direct_summation.h"
#include "verlet.h"
#include "xyz.h"
#include "thermostat.h"
#include "lattice.h"
#include "neighbors.h"
#include "lj.h"

void write_energy(std::ofstream &file, double time, double energy) {
    file << std::setw(8) << time << " " << energy << "\n";
}

int main() {
    int n = 6;
    double sigma = 2;
    double epsilon = 1;
    double m = 1;

    Atoms atoms = cubic_lattice(n, sigma);

    double end_t = 100 * std::sqrt(m * sigma * sigma / epsilon);
    double begin_t = 0;
    double step_t = 0.001 * std::sqrt(m * sigma * sigma / epsilon);
    double last_print_t = 0;
    double print_freq_t = 1 * std::sqrt(m * sigma * sigma / epsilon);
    int print_i = 0;
    std::ofstream energy_file("total_energy_001.txt");
    std::ofstream epot_file("potential_energy_001.txt");
    std::ofstream ekin_file("kinetic_energy_001.txt");

    double equi_t = step_t * 1000;
    double cutoff = sigma * 4;

    NeighborList neighbor_list;
    neighbor_list.update(atoms, cutoff);

    std::cout << "time step " << step_t << "\n";

    for (; begin_t < end_t; begin_t += step_t) {
        // compute forces
        //double e_pot = lj_direct_summation(atoms, epsilon, sigma);
        double e_pot = lj_neighbor_list(atoms, neighbor_list, cutoff, epsilon, sigma);
        double e_kin = kinetic_energy(atoms);
        double t = get_temperature(e_kin, atoms.nb_atoms());

        // apply forces
        verlet_step1(atoms, step_t, m);
        verlet_step2(atoms, step_t, m);
        double target_temp = 0.00005;
        double relaxation_t;
        if (begin_t < equi_t) {
            step_t = 0.0001 * std::sqrt(m * sigma * sigma / epsilon);
            relaxation_t = step_t * 10;
        } else {
            relaxation_t = step_t * 1000;
            step_t = 0.001 * std::sqrt(m * sigma * sigma / epsilon);
        }

        berendsen_thermostat(atoms, t, target_temp, step_t, relaxation_t);

        // compute total energy
        double e = e_pot + e_kin;
        if (begin_t - last_print_t > print_freq_t) {
            double t = get_temperature(e_kin, atoms.nb_atoms());
            std::cout << "e_pot " << std::setw(12) << e_pot
                      << " e_kin " << std::setw(12) << e_kin
                      << " t " << std::setw(4) << t << "\n";

            // log total energy
            write_energy(energy_file, begin_t, e);
            write_energy(epot_file, begin_t, e_pot);
            write_energy(ekin_file, begin_t, e_kin);

            // log positions of atoms
            last_print_t = begin_t;
            std::string num = std::to_string(print_i);
            int num_zeros = 4;
            std::string traj_file = "traj" + std::string(num_zeros - num.length(), '0') + num + ".xyz";
            print_i += 1;
            write_xyz(traj_file, atoms);

            // update neighbors
            neighbor_list.update(atoms, cutoff);
        }
    }

    energy_file.close();
    epot_file.close();
    ekin_file.close();

    return 0;
}
