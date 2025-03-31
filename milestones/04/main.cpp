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

#include <thermostat.h>

void write_energy(std::ofstream &file, double time, double energy) {
    file << std::setw(8) << time << " " << energy << "\n";
}

int main() {
    // load atoms from file
    auto [names, positions, velocities]{read_xyz_with_velocities("lj54.xyz")};
    Atoms atoms(positions);
    atoms.velocities = velocities;

    double sigma = 1;
    double epsilon = 1;
    double m = 1;

    // time step
    double begin_t = 0;
    double end_t = 100 * std::sqrt(m * sigma * sigma / epsilon);
    double step_t = 0.01 * std::sqrt(m * sigma * sigma / epsilon);
    // 0.1 bad; 0.05 bad; 0.01 good; 0.005 good;
    double print_freq_t = 1 * std::sqrt(m * sigma * sigma / epsilon);

    std::cout << "time step " << step_t << "\n";

    double last_print_t = 0;
    int print_i = 0;
    Average avg_tot; // compute averages over intervals
    std::ofstream energy_file("total_energy.txt");
    std::ofstream epot_file("potential_energy.txt");
    std::ofstream ekin_file("kinetic_energy.txt");

    for (; begin_t < end_t; begin_t += step_t) {
        // Integrator step 1
        verlet_step1(atoms, step_t, m);

        // compute forces
        double e_pot = lj_direct_summation(atoms, epsilon, sigma);

        // Integrator step 2
        verlet_step2(atoms, step_t, m);

        double e_kin = kinetic_energy(atoms);
        double e = e_pot + e_kin;
        avg_tot.add(e);

        if (begin_t - last_print_t > print_freq_t) {
            double avg_tot_r = avg_tot.result();
            double t = get_temperature_lj(e_kin, atoms.nb_atoms());
            std::cout << "e_pot " << std::setw(8) << e_pot
                      << " e_kin " << std::setw(8) << e_kin
                      << " e_tot " << std::setw(8) << avg_tot_r
                      << " t " << std::setw(4) << t << "\n";

            // log total energy
            write_energy(energy_file, begin_t, avg_tot_r);
            write_energy(epot_file, begin_t, e_pot);
            write_energy(ekin_file, begin_t, e_kin);

            // log positions of atoms
            last_print_t = begin_t;
            std::string num = std::to_string(print_i);
            int num_zeros = 4;
            std::string traj_file = "traj" + std::string(num_zeros - num.length(), '0') + num + ".xyz";
            print_i += 1;
            write_xyz(traj_file, atoms);
        }
    }

    energy_file.close();
    epot_file.close();
    ekin_file.close();

    return 0;
}
