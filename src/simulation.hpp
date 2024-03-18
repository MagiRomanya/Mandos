#ifndef SIMULABLE_H_
#define SIMULABLE_H_

#include <vector>

#include "gravity.hpp"
#include "inertia_energies.hpp"
#include "linear_algebra.hpp"
#include "particle_rigid_body_copuling.hpp"
#include "physics_state.hpp"
#include "rigid_body.hpp"
#include "fem_element.hpp"
#include "spring.hpp"
#include "particle.hpp"
#include "rod_segment.hpp"


struct Simulables {
    std::vector<Particle> particles;
    std::vector<RigidBody> rigid_bodies;
};

/** X-Macro for defining energies
 *  https://en.wikipedia.org/wiki/X_macro
 */
#define INERTIAL_ENERGY_MEMBERS \
    X(std::vector<LinearInertia>, linear_inertias) \
    X(std::vector<RotationalInertia>, rotational_inertias) \
    X(std::vector<RotationalInertiaGlobal>, rotational_global_inertias)

#define MAT(type, name) X(std::vector<FEM_Element3D<type>>, fem_elements_##name)
#define POTENTIAL_ENERGY_MEMBERS \
    X(std::vector<ParticleSpring>, particle_springs) \
    X(std::vector<RodSegment>, rod_elements) \
    FEM_MATERIAL_MEMBERS \
    X(std::vector<Gravity>, gravities)

#define X(type, name) type name;

struct Energies {
    INERTIAL_ENERGY_MEMBERS
    POTENTIAL_ENERGY_MEMBERS
};

#undef X
#undef MAT

template <typename MaterialType>
void add_FEM_element(Energies& energies, FEM_Element3D<MaterialType> element);

struct Simulation {
    Simulables simulables;
    Energies energies;

    // Integration settings
    Scalar TimeStep = 0.1;

    // Boundary conditions
    PhysicsState initial_state;
    std::vector<unsigned int> frozen_dof;
    Copulings copulings;
};

Scalar compute_energy(Scalar TimeStep, const Energies& energies, const PhysicsState& state, const PhysicsState& state0);

void compute_energy_and_derivatives(Scalar TimeStep, const Energies& energies, const PhysicsState& state, const PhysicsState& state0, EnergyAndDerivatives& out);

void simulation_step(const Simulation& simulation, PhysicsState& state);

void simulation_step(const Simulation& simulation, PhysicsState& state, EnergyAndDerivatives& out);

void update_simulation_state(const Scalar TimeStep, const Energies& energies, const Vec& dx, PhysicsState& state, const PhysicsState& state0);

#endif // SIMULABLE_H_
