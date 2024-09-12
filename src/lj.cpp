//
// Created by iter on 05.09.24.
//

#include "lj.h"

double lj_neighbor_list(Atoms &atoms, NeighborList &neighbor_list, double epsilon, double sigma) {
    // compute atoms.forces
    // return potential energy

    atoms.forces.setZero();
    double energy = 0;

    for (auto[i, k]: neighbor_list) {
        if (i < k) {
            Point_t pos_i = atoms.positions(Eigen::all, i);
            Point_t pos_k = atoms.positions(Eigen::all, k);
            Point_t r = pos_k - pos_i;
            // compute potential energy
            energy += 4 * epsilon * (std::pow(sigma, 12) * std::pow(r.norm(), -12) -
                std::pow(sigma, 6) * std::pow(r.norm(), -6));
            //std::cout << i << " " << k << " " << r.norm() << "\n";

            // compute gradient of energy - the force
            double lj = 4 * epsilon * (std::pow(sigma, 12) * (-12) * std::pow(r.norm(), -13) -
                std::pow(sigma, 6) * (-6) * std::pow(r.norm(), -7));
            Point_t f_ik = lj / r.norm() * r ; // unit vector
            atoms.forces(Eigen::all, k) -= Forces_t(f_ik);
            atoms.forces(Eigen::all, i) += Forces_t(f_ik);
        }
    }
    return energy;
}

