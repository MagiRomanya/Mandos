#ifndef PHYSICSSTATE_H_
#define PHYSICSSTATE_H_

#include <vector>
#include "linear_algebra.hpp"

struct PhysicsState {
    PhysicsState() {}
    PhysicsState(const Vec& x, const Vec& v) : x(x), v(v) {}
    Vec x;
    Vec v;

    inline void add_size(unsigned int increment_dof) {
        x.conservativeResize(x.size() + increment_dof);
        v.conservativeResize(v.size() + increment_dof);
    }

    inline unsigned int get_nDoF() const {return x.size(); }
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
