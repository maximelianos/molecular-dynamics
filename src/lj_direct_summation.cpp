//
// Copyright 2024 Maksim Velikanov
//

#include "lj_direct_summation.h"

double kinetic_energy(Atoms &atoms, double mass) {
    // e_kin = sum(m_i * v_i^2) / 2
    double energy = 0;
    for (size_t i = 0; i < atoms.nb_local; i++) {
        Point_t vel_i = atoms.velocities(Eigen::all, i);
        energy += mass * vel_i.squaredNorm();
    }
    return energy / 2;
}

double get_temperature(double e_kin, double nb_atoms) {
    // E_k = 3/2 * N_atoms k_b T -> T = 2/3 * E_k / N_atoms / k_b

    //double k_b = 1.38e+4 / 1.66 / 197; // 1.38e-23 / (1 atom mass)
    double k_b = 8.614e-5; // eV / Kelvin
    return (2.0 / 3.0) * e_kin / nb_atoms / k_b;
}

double lj_direct_summation(Atoms &atoms, double epsilon, double sigma) {
    // compute atoms.forces
    // return potential energy

    atoms.forces.setZero();
    double energy = 0;
    for (size_t i = 0; i < atoms.nb_atoms(); i++) {
        for (size_t k = i + 1; k < atoms.nb_atoms(); k++) {
            Point_t pos_i = atoms.positions(Eigen::all, i);
            Point_t pos_k = atoms.positions(Eigen::all, k);
            Point_t r = pos_k - pos_i;
            // compute potential energy
            energy += 4 * epsilon * (std::pow(sigma, 12) * std::pow(r.norm(), -12) -
                std::pow(sigma, 6) * std::pow(r.norm(), -6));
            // norm is length of vector

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
