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
#include "ducastelle.h"

void write_energy(std::ofstream &file, double time, double energy) {
    file << std::setw(8) << time << " " << energy << "\n";
}

int main() {
    double sigma = 3;
    //double epsilon = 0.01;
    double m = 197 * 103.6;

    auto [names, positions]{read_xyz("cluster_923.xyz")};
    Atoms atoms(positions);

    double begin_t = 0;
    // double end_t = 10 * std::sqrt(m * sigma * sigma / epsilon) / 10.18;
    // double step_t = 0.01 * std::sqrt(m * sigma * sigma / epsilon) / 10.18;
    double end_t = 20000;
    double step_t = 0.5;

    double last_print_t = 0;
    //double print_freq_t = 0.1 * std::sqrt(m * sigma * sigma / epsilon) / 10.18;
    double print_freq_t = 200;
    int print_i = 0;
    std::ofstream energy_file("total_energy_001.txt");
    std::ofstream epot_file("potential_energy_001.txt");
    std::ofstream ekin_file("kinetic_energy_001.txt");

    double equi_t = step_t * 100;

    NeighborList neighbor_list;
    neighbor_list.update(atoms, 10.0);

    std::cout << "time step " << step_t << "\n";

    for (; begin_t < end_t; begin_t += step_t) {
        // compute forces
        //double e_pot = lj_neighbor_list(atoms, neighbor_list, epsilon, sigma);
        //double e_pot = ducastelle(atoms, neighbor_list, sigma * 6);
        double e_pot = ducastelle(atoms, neighbor_list);
        double e_kin = kinetic_energy(atoms, m);

        // apply forces
        verlet_step1(atoms, step_t, m);
        verlet_step2(atoms, step_t, m);
        double target_temp = 500;
        double relaxation_t;
        if (begin_t < equi_t) {
            relaxation_t = step_t * 100;
        } else {
            relaxation_t = step_t * 10000;
        }

        berendsen_thermostat(atoms, target_temp, step_t, relaxation_t, m);

        // compute total energy
        double e = e_pot + e_kin;
        if (begin_t - last_print_t > print_freq_t) {
            double t = get_temperature(atoms, m);
            std::cout << "e_pot " << std::setw(12) << e_pot
                      << " e_kin " << std::setw(12) << e_kin
                      << " e_tot " << std::setw(12) << e
                      << " temp " << std::setw(4) << t << "\n";

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
            neighbor_list.update(atoms, 10.0);
        }
    }

    energy_file.close();
    epot_file.close();
    ekin_file.close();

    return 0;
}
