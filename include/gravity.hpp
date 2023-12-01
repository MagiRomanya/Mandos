#ifndef GRAVITY_H_
#define GRAVITY_H_

#include "linear_algebra.hpp"
#include "physics_state.hpp"

struct GravityParameters {
  Scalar intensity;
};

struct Gravity {
  Gravity(unsigned int index, GravityParameters params)
    : parameters(params), index(index) {}

  const GravityParameters parameters;
  const unsigned int index;

  void compute_energy_and_derivatives(const PhysicsState& state, EnergyAndDerivatives& out) const {
    const Scalar default_height = 100;
    out.energy += - parameters.intensity * (state.x(index) + default_height);
    out.gradient(index) += parameters.intensity;
    // Higher order derivatives banish
  }
};

#endif // GRAVITY_H_
