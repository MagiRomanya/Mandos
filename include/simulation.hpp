#ifndef SIMULABLE_H_
#define SIMULABLE_H_

#include <memory>
#include <vector>
#include <cassert>

#include "gravity.hpp"
#include "physics_state.hpp"
#include "rigid_body.hpp"
#include "fem_unit.hpp"
#include "spring.hpp"

struct Energies {
    std::vector<Spring> springs;
    std::vector<RigidBody> rigid_bodies;
    std::vector<FEM_Unit> fem_units;
    std::vector<Gravity> gravities;
};

struct Simulation {
    Energies energies;
    Scalar TimeStep = 0.1;
    PhysicsState initial_state;
    std::vector<Triplet> initial_mass_matrix_triplets;
    std::vector<unsigned int> frozen_dof;
};

void compute_energy_and_derivatives(const Energies& energies, const PhysicsState& state, EnergyAndDerivatives& out);

void simulation_step(const Simulation& simulation, PhysicsState& state);


#endif // SIMULABLE_H_
