//
// Created by Maksim Velikanov on 10/07/24.
//

#ifndef MD_MASTER_VERLET_H
#define MD_MASTER_VERLET_H

#include "Eigen/Core"
using Positions_t = Eigen::Array3Xd;
using Velocities_t = Eigen::Array3Xd;
using Forces_t = Eigen::Array3Xd;
using Point_t = Eigen::Vector3d;

struct Atoms {
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;

    Atoms(Positions_t &pos) :
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


void verlet_step1(Atoms &atoms, double timestep, double mass);
void verlet_step2(Atoms &atoms, double timestep, double mass);

void run_simulation(Atoms &atoms, int nb_steps);


/*void verlet_step1(double &x, double &y, double &z, double &vx, double &vy, double &vz,
                  double fx, double fy, double fz, double timestep);
void verlet_step2(double &vx, double &vy, double &vz, double fx, double fy, double fz,
                  double timestep);*/

#endif // MD_MASTER_VERLET_H
