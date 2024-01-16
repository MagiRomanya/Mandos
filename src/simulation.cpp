#include <vector>

#include "fem_unit.hpp"
#include "simulation.hpp"
#include "hard_constraints.hpp"
#include "integrators.hpp"
#include "linear_algebra.hpp"
#include "particle.hpp"
#include "physics_state.hpp"
#include "rigid_body.hpp"

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

void simulation_step(const Simulation& simulation, PhysicsState& state, EnergyAndDerivatives& out) {
    // Energy and derivatives computation
    const unsigned int nDoF = simulation.initial_state.x.size();
    EnergyAndDerivatives f(0);

#ifdef ENABLE_LAGRANGE_MULTIPLIER_CONSTRAINTS
    const unsigned int nConstraints = count_constraints(simulation.constraints);
    ConstraintsAndJacobians c(0);
#endif // ENABLE_LAGRANGE_MULTIPLIER_CONSTRAINTS
    const PhysicsState state0 = state;

    for (unsigned int i = 0; i < 1; i++) {
        f = EnergyAndDerivatives(nDoF);
#ifdef ENABLE_LAGRANGE_MULTIPLIER_CONSTRAINTS
        c = ConstraintsAndJacobians(nConstraints);
#endif // ENABLE_LAGRANGE_MULTIPLIER_CONSTRAINTS

        // Compute energy and derivatives from the energies
        // -----------------------------------------------------------------------------------------
        compute_energy_and_derivatives(simulation.TimeStep, simulation.energies, state, state0, f);

#ifdef ENABLE_LAGRANGE_MULTIPLIER_CONSTRAINTS
        // Compute constrainrs and Jacobians from the Hard Constraints
        // -----------------------------------------------------------------------------------------
        compute_constraints_and_jacobians(simulation.constraints, state, c);
#endif // ENABLE_LAGRANGE_MULTIPLIER_CONSTRAINTS

        // Integration step
        // -----------------------------------------------------------------------------------------
        Vec dx;
#ifdef ENABLE_LAGRANGE_MULTIPLIER_CONSTRAINTS
        integrate_implicit_euler(simulation, state, f, c, dx);
#else
        integrate_implicit_euler(simulation, state, f, dx);
#endif // ENABLE_LAGRANGE_MULTIPLIER_CONSTRAINTS

        // Update state
        // -----------------------------------------------------------------------------------------
        state.x_old = state0.x;
        update_simulation_state(simulation.simulables, dx, state.x);

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
