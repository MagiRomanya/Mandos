#ifndef MANDOS_GRAVITY_H_
#define MANDOS_GRAVITY_H_

#include <Mandos/linear_algebra.hpp>
#include <Mandos/physics_state.hpp>
#include <Mandos/utility_functions.hpp>

namespace mandos
{

struct GravityParameters {
    Scalar intensity;
};

struct Gravity {
    Gravity(unsigned int index, GravityParameters params)
        : parameters(params)
        , index(index)
    {
    }

    const GravityParameters parameters;
    const unsigned int index;

    inline Scalar compute_energy(Scalar TimeStep, const PhysicsState& state) const
    {
        const Scalar height = state.x(index);
        Scalar energy = -parameters.intensity * (height + default_height);
        return energy;
    }

    inline void compute_energy_gradient(Scalar TimeStep, const PhysicsState& state, Vec& grad) const
    {
        grad(index) += -parameters.intensity;  // Grad = - force
    }

    inline void compute_energy_and_derivatives(Scalar TimeStep,
                                               const PhysicsState& state,
                                               EnergyAndDerivatives& out) const
    {
        const Scalar height = state.x(index);
        out.energy += -parameters.intensity * (height + default_height);

        out.gradient(index) += -parameters.intensity;  // Grad = - force

        // ** Higher order derivatives banish **
    }

private:
    const Scalar default_height = 100;  // To avoid negative energies
};
}  // namespace mandos

#endif  // MANDOS_GRAVITY_H_
