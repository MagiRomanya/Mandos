#ifndef SPRING_H_
#define SPRING_H_

#include "linear_algebra.hpp"
#include "particle.hpp"
#include "physics_state.hpp"
#include "rigid_body.hpp"

struct SpringParameters {
    Scalar k, L0, damping;
};

struct ParticleSpring : PotentialEnergy {
    ParticleSpring(const Particle p1, const Particle p2, SpringParameters param)
        : p1(p1), p2(p2), parameters(param) {};

    const Particle p1, p2;
    const SpringParameters parameters;

    Scalar compute_energy(Scalar TimeStep, const PhysicsState& state) const;
    void compute_energy_gradient(Scalar TimeStep, const PhysicsState& state, Vec& grad) const;
    void compute_energy_and_derivatives(Scalar TimeStep, const PhysicsState& state, EnergyAndDerivatives& out) const;
};

struct RigidBodySpring : PotentialEnergy {
    const RigidBody rbA, rbB;
    const Vec3 posA, posB;
    const SpringParameters parameters;

    RigidBodySpring(const RigidBody rbA, const RigidBody rbB, Vec3 posA, Vec3 posB, SpringParameters param)
        : rbA(rbA), rbB(rbB), posA(posA), posB(posB), parameters(param) {};

    Scalar compute_energy(Scalar TimeStep, const PhysicsState& state) const;
    void compute_energy_gradient(Scalar TimeStep, const PhysicsState& state, Vec& grad) const;
    void compute_energy_and_derivatives(Scalar TimeStep, const PhysicsState& state, EnergyAndDerivatives& out) const;
};

#endif // SPRING_H_
