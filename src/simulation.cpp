#include "simulation.hpp"
#include "integrators.hpp"
#include "linear_algebra.hpp"
#include "particle.hpp"
#include "physics_state.hpp"
#include "rigid_body.hpp"
#include <vector>

void compute_energy_and_derivatives(const Energies& energies, const PhysicsState& state, EnergyAndDerivatives& out) {
    // This function is responsible of computing the energy and derivatives for each energy in the simulation.
    // Here the energies must place the energy, force and force jacobians to the correct place in the global energy and derivative structure

    // Springs
    // ---------------------------------------------------------------------
    for (size_t i = 0; i < energies.particle_springs.size(); i++) {
        energies.particle_springs[i].compute_energy_and_derivatives(state, out);
    }
    // Gravity
    // ---------------------------------------------------------------------
    for (size_t i = 0; i < energies.gravities.size(); i++) {
        energies.gravities[i].compute_energy_and_derivatives(state, out);
    }
    // FEM 3D elements
    // ---------------------------------------------------------------------
    for (size_t i = 0; i < energies.fem_elements_3d.size(); i++) {
        energies.fem_elements_3d[i].compute_energy_and_derivatives(state, out);
    }
}

void integrate(const Simulation& simulation, const PhysicsState& state, const EnergyAndDerivatives& f, Vec& new_velocities) {
    switch (simulation.integration_routine) {
    case IMPLICIT_EULER:
        integrate_implicit_euler(simulation, state, f, new_velocities);
        break;
    case SIMPLECTIC_EULER:
        integrate_simplectic_euler(simulation, state, f, new_velocities);
        break;
    }
}

void simulation_step(const Simulation& simulation, PhysicsState& state, EnergyAndDerivatives& out) {
    // Energy and derivatives computation
    const unsigned int nDoF = simulation.initial_state.x.size();
    EnergyAndDerivatives f(nDoF);

    // Compute mass matrix, kinetic energy and intrinsic energy derivatives of the simulables
    compute_simulables_energy_and_derivatives(simulation.simulables, state, f);

    // Compute energy and derivatives from the energies
    compute_energy_and_derivatives(simulation.energies, state, f);
    out = f;
    // Integration step
    Vec new_velocity;
    integrate(simulation, state, f, new_velocity);

    // Update simulables
    update_simulation_state(simulation.TimeStep, simulation.simulables, new_velocity, state);
}

void simulation_step(const Simulation& simulation, PhysicsState& state) {
    EnergyAndDerivatives f(0);
    simulation_step(simulation, state, f);
}

std::vector<Triplet> compute_global_mass_matrix(const Simulables& simulables, const PhysicsState& state) {
    std::vector<Triplet> mass_matrix_triplets;
    // PARTICLES
    // ---------------------------------------------------------------------
    for (unsigned int i = 0; i < simulables.particles.size(); i++) {
        const Particle& p = simulables.particles[i];
        // Diagonal mass matrix
        mass_matrix_triplets.emplace_back(p.index, p.index, p.mass);
        mass_matrix_triplets.emplace_back(p.index+1, p.index+1, p.mass);
        mass_matrix_triplets.emplace_back(p.index+2, p.index+2, p.mass);
    }

    // RIGID BODIES
    // ---------------------------------------------------------------------
    for (unsigned int i = 0; i < simulables.rigid_bodies.size(); i++) {
        const RigidBody& rb = simulables.rigid_bodies[i];
        // Diagonal mass matrix
        mass_matrix_triplets.emplace_back(rb.index, rb.index, rb.mass);
        mass_matrix_triplets.emplace_back(rb.index+1, rb.index+1, rb.mass);
        mass_matrix_triplets.emplace_back(rb.index+2, rb.index+2, rb.mass);

        // Rigid body moment of inertia
        const Mat3 inertia_tensor = rb.compute_inertia_tensor(state);
        for (unsigned int a = 0; a< 3; a++) {
            for (unsigned int b = 0; b< 3; b++) {
                mass_matrix_triplets.emplace_back(rb.index+3+a, rb.index+3+b, inertia_tensor(a, b));
            }
        }
    }
    return mass_matrix_triplets;
}

void compute_simulables_energy_and_derivatives(const Simulables& simulables, const PhysicsState& state, EnergyAndDerivatives& out) {
    // This function is responsible for 3 things:
    // - Computing the global mass matrix
    // - Computing the kinetic energy of the simulables
    // - Compute intrinsic forces of the simulables (the rigid body coriolis torque term)

    for (unsigned int i = 0; i < simulables.particles.size(); i++) {
        simulables.particles[i].compute_energy_and_derivatives(state, out);
    }
    for (unsigned int i = 0; i < simulables.rigid_bodies.size(); i++) {
        simulables.rigid_bodies[i].compute_energy_and_derivatives(state, out);
    }
}

void update_simulation_state(Scalar TimeStep, const Simulables& simulables, const Vec& new_velocities, PhysicsState& state) {
    for (unsigned int i = 0; i < simulables.particles.size(); i++) {
        simulables.particles[i].update_state(TimeStep, new_velocities, state);
    }
    for (unsigned int i = 0; i < simulables.rigid_bodies.size(); i++) {
        simulables.rigid_bodies[i].update_state(TimeStep, new_velocities, state);
    }
}
