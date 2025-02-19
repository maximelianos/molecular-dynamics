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
    // save energy value to a text file
    file << std::setw(8) << time << " " << energy << "\n";
}

class TextLog {
    // multiple files to save the energy

    std::ofstream energy_file_;
    std::ofstream temp_file_;
    std::ofstream stress_file_;

public:
    TextLog() {}

    // open the text files
    TextLog(std::string energy_filename, std::string temp_filename, std::string stress_filename) {
        energy_file_.open(energy_filename);
        temp_file_.open(temp_filename);
        stress_file_.open(stress_filename);
    }

    // write the values into text files
    void log(double time, double energy, double temp, double stress) {
        write_energy(energy_file_, time, energy);
        write_energy(temp_file_, time, temp);
        write_energy(stress_file_, time, stress);
    }
};


int main(int argc, char **argv) {
    // usage: ./milestone09 temp strain_diff file.xyz experiment_name
    // for example: mpirun -n 4 --oversubscribe ./milestone09 300.0 10 whisker_small.xyz 1
    MPI_Init(&argc, &argv);

    double target_temp = 300.0; // thermostat target temperature
    double strain = 5.0; // difference between beginning and end length
    std::string atoms_file = "whisker_small.xyz";
    std::string exp_name = "exp";
    if (argc >= 4) {
        target_temp = atof(argv[1]);
        strain = atof(argv[2]);
        atoms_file = argv[3];
        exp_name = argv[4];
    }

    // load gold cluster
    auto [names, positions]{read_xyz(atoms_file)}; // 923, 3871
    Atoms atoms(positions);
    int global_nb_atoms = atoms.nb_atoms();

    // automatically deduce the domain
    double dom_x = positions(0, Eigen::all).maxCoeff(); + 5.0;
    double dom_y = positions(1, Eigen::all).maxCoeff() + 5.0;

    // strain range
    double begin_strain = positions(2, Eigen::all).maxCoeff() + 1.0; // decided by EAM minimum, end must match begin!
    double end_strain = begin_strain + strain;

    // How many processors do we have?
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // initialize domain
    Domain domain(MPI_COMM_WORLD, {dom_x, dom_y, begin_strain}, {1, 1, size}, {0, 0, 1});
    int rank = domain.rank();

    // mass of a gold atom
    double m = 196.96 * 103.63; // g/mol -> [m]

    // time in femtosec
    double begin_t = 0;
    double end_t = 100000;
    double step_t = 10;

    // equilibration phase
    double equi_t = 2000;
    double print_freq_t = 1000; // 1000

    // printing status
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

    // count averages
    Average counter_temp;
    Average counter_e;
    Average counter_stress;

    double cutoff = 8.0;

    domain.enable(atoms);
    std::cout << "rank " << domain.rank() << " atoms " << atoms.nb_local << "\n";
    domain.update_ghosts(atoms, cutoff * 2);

    // update neighbours for each atom
    NeighborList neighbor_list;
    neighbor_list.update(atoms, cutoff);

    // initialize text logger
    TextLog logger;
    if (rank == 0) {
        logger = TextLog("energy.txt", "temperature.txt", "stress-" + exp_name + ".txt");
        std::cout << "time step " << step_t << "\n";
    }

    // loop over time
    for (; begin_t < end_t; begin_t += step_t) {
        // Integrator step
        verlet_step1(atoms, step_t, m);

        // compute forces between Verlet steps! Depends on ghost atoms, so we update them
        domain.exchange_atoms(atoms);
        domain.update_ghosts(atoms, cutoff * 2);
        neighbor_list.update(atoms, cutoff);

        atoms.nb_local = domain.nb_local(); // exclude ghost from e_pot
        double e_pot_local = ducastelle(atoms, neighbor_list);

        // computed stress in ducastelle
        double volume = domain.domain_length(0) * domain.domain_length(1) * domain.domain_length(2);
        double stress_local = atoms.stress(2, 2) / volume;

        // Integrator step
        verlet_step2(atoms, step_t, m);

        // Sum energies over all domain
        double e_kin_local = kinetic_energy(atoms, m); // exclude ghost from e_kin
        double e_pot = MPI::allreduce(e_pot_local, MPI_SUM, domain.communicator());
        double e_kin = MPI::allreduce(e_kin_local, MPI_SUM, domain.communicator());
        double stress = MPI::allreduce(stress_local, MPI_SUM, domain.communicator());

        // total energy
        double e = e_pot + e_kin;
        double t = get_temperature(e_kin, global_nb_atoms);

        // equilibration phase
        if (begin_t < equi_t) {
            double relaxation_t = 1000;
            berendsen_thermostat(atoms, t, target_temp, step_t, relaxation_t, m);
        }

        // compute averages
        counter_temp.add(t);
        counter_e.add(e);
        counter_stress.add(stress);

        if (begin_t > last_print_t + print_freq_t) {
            // print statistics
            last_print_t = begin_t;
            double avg_stress = counter_stress.result();
            double new_strain = (begin_t / end_t) * (end_strain - begin_strain) + begin_strain;
            if (rank == 0) {
                std::cout << "local_atoms " << std::setw(4) << domain.nb_local()
                          << " with_ghost " << std::setw(4) << atoms.nb_atoms()
                          << " e_tot " << std::setw(12) << e
                          << " temp " << std::setw(8) << t
                          << " stress " << std::setw(12) << avg_stress
                          << " new strain " << std::setw(12) << new_strain
                << "\n";
            }

            domain.disable(atoms);
            if (rank == 0) {
                logger.log(begin_t, e, t, avg_stress);
                write_xyz(get_traj_filename(), atoms);
            }
            domain.enable(atoms);

            // rescale domain, update ghost atoms

            domain.scale(atoms, {dom_x, dom_y, new_strain});
            domain.update_ghosts(atoms, cutoff * 2); // for Ducastelle
            neighbor_list.update(atoms, cutoff);
        }
    }

    MPI_Finalize();
    return 0;
}
