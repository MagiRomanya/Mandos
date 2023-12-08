#include <vector>

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

    /// Inertial energies
    // Linear
    for (size_t i = 0; i < energies.linear_inertias.size(); i++) {
        energies.linear_inertias[i].compute_energy_and_derivatives(TimeStep, state, state0, out);
    }
    // Rotation
    for (size_t i = 0; i < energies.rotational_inertias.size(); i++) {
        energies.rotational_inertias[i].compute_energy_and_derivatives(TimeStep, state, state0, out);
    }

    /// Potential energies
    // Springs
    // ---------------------------------------------------------------------
    for (size_t i = 0; i < energies.particle_springs.size(); i++) {
        energies.particle_springs[i].compute_energy_and_derivatives(TimeStep, state, out);
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

    // ... future energies here
}

void simulation_step(const Simulation& simulation, PhysicsState& state, EnergyAndDerivatives& out) {
    // Energy and derivatives computation
    const unsigned int nDoF = simulation.initial_state.x.size();
    EnergyAndDerivatives f(0);
    ConstraintsAndJacobians c;
    const PhysicsState state0 = state;

    for (unsigned int i = 0; i < 1; i++) {
        f = EnergyAndDerivatives(nDoF);
        c = ConstraintsAndJacobians();

        // Compute energy and derivatives from the energies
        // -----------------------------------------------------------------------------------------
        compute_energy_and_derivatives(simulation.TimeStep, simulation.energies, state, state0, f);

        // Compute constrainrs and Jacobians from the Hard Constraints
        // -----------------------------------------------------------------------------------------
        compute_constraints_and_jacobians(simulation.constraints, state, c);

        // Integration step
        // -----------------------------------------------------------------------------------------
        Vec dx;
        integrate_implicit_euler(simulation, state, f, c, dx);

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


void compute_constraints_and_jacobians(const HardConstraints& c, const PhysicsState& state, ConstraintsAndJacobians& out) {
    // RIGID BODY point constraints
    // ---------------------------------------------------------------------
    for (unsigned int i = 0; i < c.rb_point_constraints.size(); i++) {
        c.rb_point_constraints[i].compute_constraint_and_jacobian(state, out);
    }

    // ... future constraints here
}
