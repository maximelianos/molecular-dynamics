//
// Created by iter on 10.09.24.
//

#ifndef ATOMS_H
#define ATOMS_H

#include "Eigen/Core"

using Names_t = std::vector<std::string>;
using Positions_t = Eigen::Array3Xd;
using Velocities_t = Eigen::Array3Xd;
using Forces_t = Eigen::Array3Xd;
using Point_t = Eigen::Vector3d;

struct Atoms {
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;

    Atoms(int nb_atoms) :
          positions(3, nb_atoms),
          velocities(3, nb_atoms),
          forces(3, nb_atoms) {
        positions.setZero();
        velocities.setZero();
        forces.setZero();
    }

    Atoms(Positions_t &pos) :
          positions(pos),
          velocities(3, pos.cols()),
          forces(3, pos.cols()) {
        velocities.setZero();
        forces.setZero();
    }

    Atoms(Names_t &names, Positions_t &pos) :
          positions(pos),
          velocities(3, pos.cols()),
          forces(3, pos.cols()) {
        velocities.setZero();
        forces.setZero();
    }

    size_t nb_atoms() const {
        return positions.cols();
    }
};

#endif //ATOMS_H
