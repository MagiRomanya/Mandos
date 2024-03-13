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

    Scalar compute_energy(Scalar TimeStep, const PhysicsState& state, const PhysicsState& state0) const;
    void compute_energy_and_derivatives(Scalar TimeStep, const PhysicsState& state, const PhysicsState& state0, EnergyAndDerivatives& f) const;
    inline void update_state(const Vec& dx, Vec& x) const {
        x.segment<3>(p.index) += dx.segment<3>(p.index);
    }
};

struct RotationalInertia {
    RotationalInertia(RigidBody rb) : rb(rb) {}
    const RigidBody rb;

    Scalar compute_energy(Scalar TimeStep, const PhysicsState& state, const PhysicsState& state0) const;
    void compute_energy_and_derivatives(Scalar TimeStep, const PhysicsState& state, const PhysicsState& state0, EnergyAndDerivatives& f) const;
    void update_state(const Vec& dx, Vec& x) const;
};

struct RotationalInertiaGlobal {
    RotationalInertiaGlobal(RigidBody rb) : rb(rb) {}
    const RigidBody rb;

    Scalar compute_energy(Scalar TimeStep, const PhysicsState& state, const PhysicsState& state0) const;
    void compute_energy_and_derivatives(Scalar TimeStep, const PhysicsState& state, const PhysicsState& state0, EnergyAndDerivatives& f) const;
    void update_state(const Vec& dx, Vec& x) const;
};



struct Simulation;
void add_particle_to_simulation(Simulation& simulation, const Particle& p);
void add_rigid_body_to_simulation(Simulation& simulation, const RigidBody& rb, const bool global = false);

#endif // INERTIA_ENERGIES_H_
