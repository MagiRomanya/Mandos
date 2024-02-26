#ifndef GRAVITY_H_
#define GRAVITY_H_

#include "linear_algebra.hpp"
#include "physics_state.hpp"
#include "utility_functions.hpp"

struct GravityParameters {
  Scalar intensity;
};

struct Gravity {
  Gravity(unsigned int index, GravityParameters params)
    : parameters(params), index(index) {}

  const GravityParameters parameters;
  const unsigned int index;

  inline Scalar compute_energy(Scalar TimeStep, const PhysicsState& state) const {
    const Scalar height = state.x(index);
    Scalar energy = - parameters.intensity * (height + default_height);
    return energy;
  }

  inline void compute_energy_and_derivatives(Scalar TimeStep, const PhysicsState& state, EnergyAndDerivatives& out) const {
    const Scalar height = state.x(index);
    out.energy += - parameters.intensity * (height + default_height);

    out.gradient(index) += - parameters.intensity; // Grad = - force

    // ** Higher order derivatives banish **
  }
  private:
    const Scalar default_height = 100; // To avoid negative energies
};

#endif // GRAVITY_H_
