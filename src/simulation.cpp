#include <vector>

#include "fem_unit.hpp"
#include "simulation.hpp"
#include "integrators.hpp"
#include "linear_algebra.hpp"
#include "particle.hpp"
#include "physics_state.hpp"
#include "rigid_body.hpp"
#include "utility_functions.hpp"

void compute_energy_and_derivatives(Scalar TimeStep, const Energies& energies, const PhysicsState& state, const PhysicsState& state0, EnergyAndDerivatives& out) {
    // This function is responsible of computing the energy and derivatives for each energy in the simulation.
    // Here the energies must place the energy, gradient and Hessians to the correct place in the global energy and derivative structure

#define X(type, energy) \
    for (size_t i = 0; i < energies.energy.size(); i++) { \
        energies.energy[i].compute_energy_and_derivatives(TimeStep, state, state0, out); \
    }
    INERTIAL_ENERGY_MEMBERS
#undef X

#define MAT(type, name) X(std::vector<FEM_Element3D<type>>, fem_elements_##name)
#define X(type, energy) \
    for (size_t i = 0; i < energies.energy.size(); i++) { \
        energies.energy[i].compute_energy_and_derivatives(TimeStep, state, out); \
    }
    POTENTIAL_ENERGY_MEMBERS
#undef X
#undef MAT
}

Scalar compute_energy(Scalar TimeStep, const Energies& energies, const PhysicsState& state, const PhysicsState& state0) {
    // This function is responsible of computing the energy and derivatives for each energy in the simulation.
    // Here the energies must place the energy, gradient and Hessians to the correct place in the global energy and derivative structure
    Scalar phi = 0;
#define X(type, energy) \
    for (size_t i = 0; i < energies.energy.size(); i++) { \
        phi += energies.energy[i].compute_energy(TimeStep, state, state0); \
    }
    INERTIAL_ENERGY_MEMBERS
#undef X

#define MAT(type, name) X(std::vector<FEM_Element3D<type>>, fem_elements_##name)
#define X(type, energy) \
    for (size_t i = 0; i < energies.energy.size(); i++) { \
        phi += energies.energy[i].compute_energy(TimeStep, state); \
    }
    POTENTIAL_ENERGY_MEMBERS
#undef X
#undef MAT
        return phi;
}

void simulation_step(const Simulation& simulation, PhysicsState& state, EnergyAndDerivatives& out) {
    // Energy and derivatives computation
    const unsigned int nDoF = simulation.initial_state.x.size();
    EnergyAndDerivatives f(0);

    const PhysicsState state0 = state;
    // const Vec dx0 = (state.x - state.x_old); // h * v
    // update_simulation_state(simulation.simulables, dx0, state.x);
    // update_simulation_state(simulation.simulables, dx0, state.x_old);

    for (unsigned int i = 0; i < 1; i++) {

        // Compute energy and derivatives from the energies
        // -----------------------------------------------------------------------------------------
        f = EnergyAndDerivatives(nDoF);
        compute_energy_and_derivatives(simulation.TimeStep, simulation.energies, state, state0, f);
        Scalar energy = f.energy;

        // Integration step
        // -----------------------------------------------------------------------------------------
        Vec dx;
        integrate_implicit_euler(simulation, state, f, dx);

        // Update state
        // -----------------------------------------------------------------------------------------
        const PhysicsState stepState = state;
        state.x_old = state0.x;
        update_simulation_state(simulation.simulables, dx, state.x);

        Scalar lineSearchEnergy = compute_energy(simulation.TimeStep, simulation.energies, state, state0);
        Scalar alpha = 1.0f;
        const bool enableLineSearch = false;
        while (lineSearchEnergy > energy && enableLineSearch) {
            DEBUG_LOG(alpha);
            alpha /= 2.0f;
            state = stepState;
            update_simulation_state(simulation.simulables, alpha * dx, state.x);
            lineSearchEnergy = compute_energy(simulation.TimeStep, simulation.energies, state, state0);
        }

    }
    // Output energy derivatives
    out = f;
}

void simulation_step(const Simulation& simulation, PhysicsState& state) {
    EnergyAndDerivatives f(0);
    simulation_step(simulation, state, f);
}

void update_simulation_state(const Simulables& simulables, const Vec& dx, Vec& x) {
    for (unsigned int i = 0; i < simulables.particles.size(); i++) {
        simulables.particles[i].update_state(dx, x);
    }
    for (unsigned int i = 0; i < simulables.rigid_bodies.size(); i++) {
        simulables.rigid_bodies[i].update_state(dx, x);
    }
}

#define MAT(type, name) template void add_FEM_element(Energies& energies, FEM_Element3D<type> element);
FEM_MATERIAL_MEMBERS
#undef MAT

template <typename MaterialType>
void add_FEM_element(Energies& energies, FEM_Element3D<MaterialType> element) {
#define MAT(type, name)                                     \
        if constexpr (std::is_same<MaterialType, type>()) { \
            energies.fem_elements_##name.push_back(element); \
        }

    FEM_MATERIAL_MEMBERS

#undef MAT
}
