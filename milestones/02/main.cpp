#include "verlet.h"
#include <Eigen/Dense>
#include <iostream>
#include <filesystem>

#ifdef USE_MPI
#include <mpi.h>
#endif


int main(int argc, char *argv[]) {
    int rank = 0, size = 1;

    // Below is some MPI code, try compiling with `cmake -DUSE_MPI=ON ..`
#ifdef USE_MPI
    MPI_Init(&argc, &argv);

    // Retrieve process infos
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

    std::cout << "Hello from milestone 2\n";
    std::cout << "Hello I am rank " << rank << " of " << size << "\n";

    if (rank == 0) {
        Positions_t pos(3, 1);
        pos(0, 0) = 1;
        pos(1, 0) = 0;
        pos(2, 0) = 0;
        Atoms atoms(pos);
        atoms.forces(0, 0) = 1;
        atoms.forces(1, 0) = 0;
        atoms.forces(2, 0) = 0;

        std::cout << atoms.positions << "\n";

        int nb_steps = 10;
        for (int i = 0; i < nb_steps; i++) {
            std::cout << "Step: " << i << "\n";
            verlet_step1(atoms, 1, 1);
            verlet_step2(atoms, 1, 1);
            std::cout << atoms.positions << "\n";
        }
    }

#ifdef USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
