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
#include "mpi.h"

#include "lj_direct_summation.h"
#include "verlet.h"
#include "xyz.h"
#include "thermostat.h"
#include "lattice.h"
#include "neighbors.h"
#include "lj.h"
#include "ducastelle.h"

#include "mpi_support.h"
#include "domain.h"

void write_energy(std::ofstream &file, double time, double energy) {
    file << std::setw(8) << time << " " << energy << "\n";
}

class TextLog {
    std::ofstream energy_file_;
    std::ofstream temp_file_;

public:
    TextLog() {}

    TextLog(char *energy_filename, char *temp_filename) {
        energy_file_.open(energy_filename);
        temp_file_.open(temp_filename);
    }

    void log(double time, double energy, double temp) {
        write_energy(energy_file_, time, energy);
        write_energy(temp_file_, time, temp);
    }
};


int main(int argc, char **argv) {
    // usage: ./milestone08 lx ly lz nx ny nz
    MPI_Init(&argc, &argv);

    float lx = 60.0;
    float ly = 60.0;
    float lz = 60.0;
    int nx = 1;
    int ny = 1;
    int nz = 1;
    if (argc >= 7) {
        lx = atof(argv[1]);
        ly = atof(argv[2]);
        lz = atof(argv[3]);
        nx = atoi(argv[4]);
        ny = atoi(argv[5]);
        nz = atoi(argv[6]);
    }

    Domain domain(MPI_COMM_WORLD, {lx, ly, lz}, {nx, ny, nz}, {0, 0, 0});
    int rank = domain.rank();

    double m = 196.96 * 103.63; // g/mol -> [m]

    // load gold cluster
    auto [names, positions]{read_xyz("cluster_3871.xyz")}; // 923, 3871
    Atoms atoms(positions);
    for (size_t i = 0; i < atoms.nb_atoms(); i++) {
        Point_t pos_i = atoms.positions(Eigen::all, i);
        // add offset from 0
        atoms.positions(Eigen::all, i) = pos_i + Eigen::Vector3d(5.0, 5.0, 5.0);
    }
    int global_nb_atoms = atoms.nb_atoms();

    // time in femtosec
    double begin_t = 0;
    double end_t = 10000;
    double step_t = 10;

    // equilibration phase
    double equi_t = 500;
    double print_freq_t = 100; // 1000

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

    double temp_sum = 0;
    double e_tot_sum = 0;

    double cutoff = 9.0;

    domain.enable(atoms);
    std::cout << "rank " << domain.rank() << " atoms " << atoms.nb_local << "\n";
    domain.update_ghosts(atoms, cutoff * 2);
    NeighborList neighbor_list;
    neighbor_list.update(atoms, cutoff);

    TextLog logger;
    if (rank == 0) {
        logger = TextLog("energy.txt", "temperature.txt");
        std::cout << "time step " << step_t << "\n";
    }

    for (; begin_t < end_t; begin_t += step_t) {
        // Integrator step, independant
        verlet_step1(atoms, step_t, m);

        // compute forces between Verlet steps! Depends on ghost atoms
        domain.exchange_atoms(atoms);
        domain.update_ghosts(atoms, cutoff * 2);
        neighbor_list.update(atoms, cutoff);

        atoms.nb_local = domain.nb_local(); // exclude ghost from e_pot
        double e_pot_local = ducastelle(atoms, neighbor_list);

        // Integrator step, independant
        verlet_step2(atoms, step_t, m);

        // Sum energies over all domain
        double e_kin_local = kinetic_energy(atoms, m); // exclude ghost from e_kin
        double e_pot = MPI::allreduce(e_pot_local, MPI_SUM, domain.communicator());
        double e_kin = MPI::allreduce(e_kin_local, MPI_SUM, domain.communicator());

        double e = e_pot + e_kin;
        double t = get_temperature(e_kin, global_nb_atoms);

        // berendsen thermostat not needed
        if (begin_t < equi_t) {
            double target_temp = 300;
            double relaxation_t = 1000;
            //berendsen_thermostat(atoms, target_temp, step_t, relaxation_t, m);
        }

        temp_sum += t;
        e_tot_sum += e;

        if (begin_t > last_print_t + print_freq_t) {
            last_print_t = begin_t;

            if (rank == 0) {
                std::cout << "local_atoms " << std::setw(4) << domain.nb_local()
                          << " with_ghost " << std::setw(4) << atoms.nb_atoms()
                          << " e_pot " << std::setw(12) << e_pot
                          << " e_kin " << std::setw(12) << e_kin
                          << " e_tot " << std::setw(12) << e
                          << " temp " << std::setw(4) << t << "\n";
            }

            domain.disable(atoms);
            if (rank == 0) {
                logger.log(begin_t, e, t);
                write_xyz(get_traj_filename(), atoms);
            }
            domain.enable(atoms);
            domain.update_ghosts(atoms, cutoff * 2); // for Ducastelle
            neighbor_list.update(atoms, cutoff);
        }
    }

    MPI_Finalize();
    return 0;
}
