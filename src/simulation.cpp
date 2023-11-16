#include "simulation.hpp"
#include "integrators.hpp"
#include "linear_algebra.hpp"
#include "particle.hpp"
#include "physics_state.hpp"
#include "rigid_body.hpp"
#include <vector>

void compute_energy_and_derivatives(const Energies& energies, const PhysicsState& state, EnergyAndDerivatives& out) {

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
}

void simulation_step(const Simulation& simulation, PhysicsState& state) {
    // Energy and derivatives computation
    const unsigned int nDoF = simulation.initial_state.x.size();
    EnergyAndDerivatives f(nDoF);
    compute_energy_and_derivatives(simulation.energies, state, f);

    // Integration step
    integrate_implicit_euler(simulation, &state, f);
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
                mass_matrix_triplets.emplace_back(rb.index+a, rb.index+b, inertia_tensor(a, b));
            }
        }
    }
    return mass_matrix_triplets;
}
