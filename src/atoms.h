/*
 * Copyright 2024 Maxim Velikanov
 * MIT license
 *
 * The Atoms class holds the atomic system positions, velocities and forces.
 */

#ifndef ATOMS_H
#define ATOMS_H

#include "Eigen/Core"

using Names_t = std::vector<std::string>;
using Positions_t = Eigen::Array3Xd;
using Velocities_t = Eigen::Array3Xd;
using Forces_t = Eigen::Array3Xd;
using Masses_t = Eigen::ArrayXd;
using Point_t = Eigen::Vector3d;

struct Atoms {
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;
    Masses_t masses;
    int nb_local; // domain decomposition
    Eigen::Matrix3d stress; // whisker stretching

    Atoms(int nb_atoms) :
          positions(3, nb_atoms),
          velocities(3, nb_atoms),
          forces(3, nb_atoms),
          masses(nb_atoms) {
        positions.setZero();
        velocities.setZero();
        forces.setZero();
        masses.setZero();
        nb_local = nb_atoms;
    }

    Atoms(Positions_t &pos) :
          positions(pos),
          velocities(3, pos.cols()),
          forces(3, pos.cols()),
          masses(pos.cols()) {
        velocities.setZero();
        forces.setZero();
        masses.setZero();
        nb_local = pos.cols();
    }

    Atoms(Names_t &names, Positions_t &pos) :
          positions(pos),
          velocities(3, pos.cols()),
          forces(3, pos.cols()),
          masses(pos.cols()) {
        velocities.setZero();
        forces.setZero();
        masses.setZero();
        nb_local = pos.cols();
    }

    size_t nb_atoms() const {
        return positions.cols();
    }

    // Domain decomposition. Change the size of arrays for ghost atoms
    // or to distribute the subdomains.
    void resize(size_t n) {
        positions.conservativeResize(3, n);
        velocities.conservativeResize(3, n);
        forces.conservativeResize(3, n);
        masses.conservativeResize(n);
        nb_local = n;
    };
};

#endif //ATOMS_H
