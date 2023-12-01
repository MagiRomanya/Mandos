#ifndef INERTIA_ENERGIES_H_
#define INERTIA_ENERGIES_H_


#include "linear_algebra.hpp"
#include "particle.hpp"
#include "physics_state.hpp"
#include "rigid_body.hpp"

struct LinearInertia {
    LinearInertia(Particle p) : Mass(Mat3::Identity() * p.mass), p(p) {}

    const Mat3 Mass;
    const Particle p;

    void compute_energy_and_derivatives(Scalar TimeStep, const PhysicsState& state, EnergyAndDerivatives& f) const;
};

struct RotationalInertia {
    RotationalInertia(RigidBody rb) : rb(rb) {}
    const RigidBody rb;
    void compute_energy_and_derivatives(Scalar TimeStep, const PhysicsState& state, EnergyAndDerivatives& f) const;
};


#endif // INERTIA_ENERGIES_H_
