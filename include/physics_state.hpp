#ifndef PHYSICSSTATE_H_
#define PHYSICSSTATE_H_

#include <vector>

#include "linear_algebra.hpp"

struct PhysicsState {
    Vec x;
    Vec v; // (x - x_old) / TimeStep

    void add_size(unsigned int increment_dof) {
        x.conservativeResize(x.size() + increment_dof);
        v.conservativeResize(v.size() + increment_dof);
    }
};

struct EnergyAndDerivatives {
    EnergyAndDerivatives(unsigned int nDoF) {
        energy = 0;
        force.setZero(nDoF);
        df_dx_triplets.clear();
    }
    // Container
    Scalar energy;
    Vec force;
    std::vector<Triplet> df_dx_triplets;
};


#endif // PHYSICSSTATE_H_
