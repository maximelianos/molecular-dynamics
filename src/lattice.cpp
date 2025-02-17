//
// Copyright 2024 Maksim Velikanov
//

#include "lattice.h"

Atoms cubic_lattice(int n, double sigma) {
    Atoms atoms(n*n*n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                atoms.positions(0, i*n*n + j*n + k) = i * sigma * 1.13;
                atoms.positions(1, i*n*n + j*n + k) = j * sigma * 1.13;
                atoms.positions(2, i*n*n + j*n + k) = k * sigma * 1.13;
            }
        }
    }
    return atoms;
}