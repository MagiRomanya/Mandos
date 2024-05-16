#ifndef MANDOS_PHYSICSSTATE_H_
#define MANDOS_PHYSICSSTATE_H_

#include <vector>

#include <Mandos/linear_algebra.hpp>

namespace mandos
{
struct PhysicsState {
    PhysicsState()
    {
    }
    PhysicsState(const Vec& x, const Vec& v)
        : x(x)
        , v(v)
    {
    }
    Vec x;
    Vec v;

    inline void add_size(unsigned int increment_dof)
    {
        x.conservativeResize(x.size() + increment_dof);
        v.conservativeResize(v.size() + increment_dof);
    }

    inline unsigned int get_nDoF() const
    {
        return static_cast<unsigned int>(x.size());
    }
};

struct EnergyAndDerivatives {
    EnergyAndDerivatives()
    {
    }
    EnergyAndDerivatives(unsigned int nDoF)
    {
        energy = 0;
        gradient.setZero(nDoF);
        hessian_triplets.clear();
    }
    inline void clear(unsigned int nDoF)
    {
        energy = 0;
        gradient.setZero(nDoF);
        hessian_triplets.clear();
    }
    // Container
    Scalar energy;
    Vec gradient;
    std::vector<Triplet> hessian_triplets;
};

}  // namespace mandos

#endif  // MANDOS_PHYSICSSTATE_H_
