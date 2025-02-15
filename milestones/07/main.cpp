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


void run_heat_capacity() {
    double m = 196.96 * 103.63; // g/mol -> [m]

    // load gold cluster
    auto [names, positions]{read_xyz("whisker.xyz")}; // 923, 3871
    Atoms atoms(positions);

    // Small cluster
    // end t = 100 000
    // step_t = 10
    // deposit_t = 4000
    // delta_q = 20

    // Big cluster
    // delta_q = 100

    // time in femtosec
    double begin_t = 0;
    double end_t = 100000;
    double step_t = 10;

    // equilibration phase
    double equi_t = 500;

    double heat_deposit_t = 4000; // 1000
    double print_freq_t = 1000; // 100

    double last_print_t = 0;
    int print_i = 0;
    // log positions of atoms
    auto get_traj_filename = [&print_i] -> std::string {
        std::string num = std::to_string(print_i);
        const int num_zeros = 4;
        std::string traj_file = "traj" + std::string(num_zeros - num.length(), '0') + num + ".xyz";
        print_i += 1;
        return traj_file;
    };

    std::ofstream energy_file("total_energy.txt");
    std::ofstream epot_file("potential_energy.txt");
    std::ofstream ekin_file("kinetic_energy.txt");
    std::ofstream temp_file("temperature.txt");

    std::ofstream interval_energy_file("interval_e_total.txt");
    std::ofstream interval_temp_file("interval_temperature.txt");

    double last_heat_t = 0;
    double temp_sum = 0;
    double e_tot_sum = 0;

    double cutoff = 10.0;
    NeighborList neighbor_list;
    neighbor_list.update(atoms, cutoff);

    std::cout << "time step " << step_t << "\n";

    for (; begin_t < end_t; begin_t += step_t) {
        verlet_step1(atoms, step_t, m);
        // compute forces between Verlet steps!
        double e_pot = ducastelle(atoms, neighbor_list);
        verlet_step2(atoms, step_t, m);

        double e_kin = kinetic_energy(atoms, m);
        double e = e_pot + e_kin;
        double t = get_temperature(atoms);

        double target_temp = 300;
        double relaxation_t;
        if (begin_t < equi_t) {
            relaxation_t = 1000;
            berendsen_thermostat(atoms, target_temp, step_t, relaxation_t, m);
        }

        // deposit heat
        temp_sum += t;
        e_tot_sum += e;
        if (last_heat_t + heat_deposit_t < begin_t) {
            double steps = heat_deposit_t / step_t;
            std::cout << "DEPOSIT"
                      << " e_avg " << std::setw(12) << e_tot_sum / steps
                      << " t_avg " << std::setw(12) << temp_sum / steps
                    << "\n";

            double t_avg = temp_sum / steps;
            double e_avg = e_tot_sum / steps;
            write_energy(interval_energy_file, begin_t, e_avg);
            write_energy(interval_temp_file, begin_t, t_avg);
            write_xyz(get_traj_filename(), atoms);

            double delta_q = 200.0;
            double lambda = std::sqrt(1 + delta_q / e_kin); // add constant heat
            atoms.velocities = atoms.velocities * lambda;

            last_heat_t = begin_t;
            temp_sum = 0;
            e_tot_sum = 0;
        }


        if (begin_t - last_print_t > print_freq_t) {
            std::cout << "e_pot " << std::setw(12) << e_pot
                      << " e_kin " << std::setw(12) << e_kin
                      << " e_tot " << std::setw(12) << e
                      << " temp " << std::setw(4) << t << "\n";
            last_print_t = begin_t;

            // log total energy
            write_energy(energy_file, begin_t, e);
            write_energy(epot_file, begin_t, e_pot);
            write_energy(ekin_file, begin_t, e_kin);
            write_energy(temp_file, begin_t, t);
            //write_xyz(get_traj_filename(), atoms);

            // update neighbors
            neighbor_list.update(atoms, cutoff);
        }
    }

    energy_file.close();
    epot_file.close();
    ekin_file.close();
    interval_energy_file.close();
    interval_temp_file.close();
}

int main() {
    run_heat_capacity();
    return 0;
}
