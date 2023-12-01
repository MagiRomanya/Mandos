#ifndef PHYSICSSTATE_H_
#define PHYSICSSTATE_H_

#include <vector>

#include "linear_algebra.hpp"

struct PhysicsState {
    Vec x;
    Vec v; // v = x - x_old
    Vec x_old;
    Vec v_old;

    void add_size(unsigned int increment_dof) {
        x.conservativeResize(x.size() + increment_dof);
        v.conservativeResize(v.size() + increment_dof);
    }
};

struct EnergyAndDerivatives {
    EnergyAndDerivatives(unsigned int nDoF) {
        energy = 0;
        jacobian.setZero(nDoF);
        hessian_triplets.clear();
    }
    // Container
    Scalar energy;
    Vec jacobian;
    std::vector<Triplet> hessian_triplets;
};


#endif // PHYSICSSTATE_H_
