#ifndef SIMULABLE_H_
#define SIMULABLE_H_

#include <memory>
#include <vector>
#include <cassert>

#include "gravity.hpp"
#include "linear_algebra.hpp"
#include "physics_state.hpp"
#include "rigid_body.hpp"
#include "fem_unit.hpp"
#include "spring.hpp"
#include "particle.hpp"

struct Simulables {
    std::vector<Particle> particles;
    std::vector<RigidBody> rigid_bodies;
};

struct Energies {
    std::vector<ParticleSpring> particle_springs;
    std::vector<FEM_Unit> fem_units;
    std::vector<Gravity> gravities;
};

struct Simulation {
    Simulables simulables;
    Energies energies;
    Scalar TimeStep = 0.1;
    PhysicsState initial_state;
    std::vector<unsigned int> frozen_dof;
};

void compute_energy_and_derivatives(const Energies& energies, const PhysicsState& state, EnergyAndDerivatives& out);

std::vector<Triplet> compute_global_mass_matrix(const Simulables& simulables, const PhysicsState& state);

void simulation_step(const Simulation& simulation, PhysicsState& state);


#endif // SIMULABLE_H_
