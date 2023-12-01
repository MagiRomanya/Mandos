#ifndef PHYSICSSTATE_H_
#define PHYSICSSTATE_H_

#include <vector>

#include "linear_algebra.hpp"

struct PhysicsState {
    Vec x;
    Vec x_old;
    Vec x_old2;

    void add_size(unsigned int increment_dof) {
        x.conservativeResize(x.size() + increment_dof);
        x_old.conservativeResize(x_old.size() + increment_dof);
        x_old2.conservativeResize(x_old2.size() + increment_dof);
    }
};

struct EnergyAndDerivatives {
    EnergyAndDerivatives(unsigned int nDoF) {
        energy = 0;
        gradient.setZero(nDoF);
        hessian_triplets.clear();
    }
    // Container
    Scalar energy;
    Vec gradient;
    std::vector<Triplet> hessian_triplets;
};


#endif // PHYSICSSTATE_H_
