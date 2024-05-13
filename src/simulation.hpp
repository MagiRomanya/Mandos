#ifndef SIMULABLE_H_
#define SIMULABLE_H_

#include <vector>

#include "colliders.hpp"
#include "linear_algebra.hpp"
#include "particle_rigid_body_coupling.hpp"
#include "physics_state.hpp"
#include "rigid_body.hpp"
#include "particle.hpp"
#include "energies.hpp"


struct Simulables {
    std::vector<Particle> particles;
    std::vector<RigidBody> rigid_bodies;
};

struct Simulation {
    Simulables simulables;
    Energies energies;
    Colliders colliders;

    // Integration settings
    Scalar TimeStep = 0.1;
    unsigned int MaxNewtonIterations = 1;
    bool enable_line_search = false;

    // Boundary conditions
    PhysicsState initial_state;
    std::vector<unsigned int> frozen_dof;
    Couplings couplings;
};

Scalar compute_energy(Scalar TimeStep, const Energies& energies, const PhysicsState& state, const PhysicsState& state0);

Vec compute_energy_gradient(Scalar TimeStep, const Energies& energies, const PhysicsState& state, const PhysicsState& state0);

void compute_energy_and_derivatives(Scalar TimeStep, const Energies& energies, const PhysicsState& state, const PhysicsState& state0, EnergyAndDerivatives& out);

void simulation_step(const Simulation& simulation, PhysicsState& state);

void simulation_step(const Simulation& simulation, PhysicsState& state, EnergyAndDerivatives& out);

void update_simulation_state(const Scalar TimeStep, const Energies& energies, const Vec& dx, PhysicsState& state, const PhysicsState& state0);

#endif // SIMULABLE_H_
