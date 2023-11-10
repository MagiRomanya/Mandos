#ifndef PHYSICSSTATE_H_
#define PHYSICSSTATE_H_

#include <vector>

#include "linear_algebra.hpp"

struct PhysicsState {
    Vec x;
    Vec v; // (x - x_old) / TimeStep

    void resize(unsigned int increment_dof) {
        x.conservativeResize(x.size() + increment_dof);
        v.conservativeResize(v.size() + increment_dof);
    }
};

struct EnergyAndDerivatives {
    // Container
    Scalar energy;
    Vec force;
    SparseMat df_dx;
    std::vector<Triplet> df_dx_triplets;
};


#endif // PHYSICSSTATE_H_
