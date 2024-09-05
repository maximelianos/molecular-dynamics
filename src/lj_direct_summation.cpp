//
// Created by iter on 05.09.24.
//

#include <iostream>
#include "verlet.h"
#include "lj_direct_summation.h"

double lj_direct_summation(Atoms &atoms, double epsilon, double sigma) {
    atoms.forces.setZero();
    for (size_t i = 0; i < atoms.nb_atoms(); i++) {
        for (size_t j = i + 1; j < atoms.nb_atoms(); j++) {
            Point_t pos_i = atoms.positions(Eigen::all, i);
            Point_t pos_j = atoms.positions(Eigen::all, j);
            Point_t r = pos_j - pos_i;
            //std::cout << i << " " << j << " " << r.norm() << "\n";

            float lj = 4 * epsilon * (std::pow(sigma, 12) * (-12) * std::pow(r.norm(), -13) -
                std::pow(sigma, 6) * (-6) * std::pow(r.norm(), -7));
            Point_t f_ik = lj / r.norm() * r ; // unit vector
            atoms.forces(Eigen::all, j) += Forces_t(f_ik);
            atoms.forces(Eigen::all, i) -= (Forces_t) f_ik;
        }
    }
    return 0;
}