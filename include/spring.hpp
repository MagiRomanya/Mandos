#ifndef SPRING_H_
#define SPRING_H_

#include "linear_algebra.hpp"
#include "particle.hpp"
#include "physics_state.hpp"

struct SpringParameters {
    Scalar k, L0, damping;

    Scalar get_energy(Scalar L) const;
    Vec3 get_energy_gradient(const Vec3& x1, const Vec3& x2, const Vec3& v1, const Vec3& v2, Scalar L) const;
    Mat3 get_energy_hessian(Scalar TimeStep, const Vec3& x1, const Vec3& x2, Scalar L) const;
};

struct ParticleSpring {
    ParticleSpring(const Particle p1, const Particle p2, SpringParameters param)
        : p1(p1), p2(p2), parameters(param) {};

    const Particle p1;
    const Particle p2;
    const SpringParameters parameters;

    Scalar compute_energy(Scalar TimeStep, const PhysicsState& state) const;
    void compute_energy_and_derivatives(Scalar TimeStep, const PhysicsState& state, EnergyAndDerivatives& out) const;
};

#endif // SPRING_H_
