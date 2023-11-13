#include "simulation.hpp"
#include "integrators.hpp"
#include "physics_state.hpp"

void compute_energy_and_derivatives(const Energies& energies, const PhysicsState& state, EnergyAndDerivatives& out) {
    for (size_t i = 0; i < energies.springs.size(); i++) {
        energies.springs[i].compute_energy_and_derivatives(state, out);
    }

    for (size_t i = 0; i < energies.gravities.size(); i++) {
        energies.gravities[i].compute_energy_and_derivatives(state, out);
    }

    // for (size_t i = 0; i < energies.rigid_bodies.size(); i++) {
    //     energies.rigid_bodies[i].compute_energy_and_derivatives(state, out);
    // }
    //
    // ...
}

void simulation_step(const Simulation& simulation, PhysicsState& state) {
    // Energy and derivatives computation
    const unsigned int nDoF = simulation.initial_state.x.size();
    EnergyAndDerivatives f(nDoF);
    compute_energy_and_derivatives(simulation.energies, state, f);

    // Integration step
    integrate_implicit_euler(simulation, &state, f);
}
