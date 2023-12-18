#ifndef SIMULABLE_H_
#define SIMULABLE_H_

#include <memory>
#include <vector>
#include <cassert>

#include "gravity.hpp"
#include "hard_constraints.hpp"
#include "inertia_energies.hpp"
#include "linear_algebra.hpp"
#include "particle_rigid_body_copuling.hpp"
#include "physics_state.hpp"
#include "rigid_body.hpp"
#include "fem_unit.hpp"
#include "spring.hpp"
#include "particle.hpp"


struct Simulables {
    std::vector<Particle> particles;
    std::vector<RigidBody> rigid_bodies;
};

/* X-Macro for defining energies */
#define INERTIAL_ENERGY_MEMBERS \
    X(std::vector<LinearInertia>, linear_inertias) \
    X(std::vector<RotationalInertia>, rotational_inertias)

#define POTENTIAL_ENERGY_MEMBERS \
    X(std::vector<ParticleSpring>, particle_springs) \
    X(std::vector<FEM_Element3D>, fem_elements_3d) \
    X(std::vector<Gravity>, gravities)

struct Energies {
#define X(type, name) type name;
    INERTIAL_ENERGY_MEMBERS
    POTENTIAL_ENERGY_MEMBERS
#undef X
};

struct Simulation {
    Simulables simulables;
    Energies energies;

    // Integration settings
    Scalar TimeStep = 0.01;

    // Boundary conditions
    PhysicsState initial_state;
    std::vector<unsigned int> frozen_dof;
    std::vector<ParticleRigidBodyCopuling> copulings;

#ifdef ENABLE_LAGRANGE_MULTIPLIER_CONSTRAINTS
    HardConstraints constraints;
#endif
};

void compute_energy_and_derivatives(Scalar TimeStep, const Energies& energies, const PhysicsState& state, const PhysicsState& state0, EnergyAndDerivatives& out);

void simulation_step(const Simulation& simulation, PhysicsState& state);

void simulation_step(const Simulation& simulation, PhysicsState& state, EnergyAndDerivatives& out);

void update_simulation_state(const Simulables& simulables, const Vec& dx, Vec& x);

#endif // SIMULABLE_H_
