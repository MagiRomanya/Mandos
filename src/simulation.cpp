#include <cmath>
#include <vector>

#include "clock.hpp"
#include "colliders.hpp"
#include "simulation.hpp"
#include "integrators.hpp"
#include "linear_algebra.hpp"
#include "physics_state.hpp"
#include "utility_functions.hpp"

void compute_energy_and_derivatives_finite(Scalar TimeStep, const Energies& energies, const PhysicsState& state, const PhysicsState& state0, EnergyAndDerivatives& out);

void compute_energy_and_derivatives(Scalar TimeStep, const Energies& energies, const PhysicsState& state, const PhysicsState& state0, EnergyAndDerivatives& out) {
    // This function is responsible of computing the energy and derivatives for each energy in the simulation.
    // Here the energies must place the energy, gradient and Hessians to the correct place in the global energy and derivative structure
    energies.inertial_energies.for_each(ComputeInertialEnergyAndDerivativesVisitor(TimeStep, state, state0, out));
    energies.potential_energies.for_each(ComputePotentialEnergyAndDerivativesVisitor(TimeStep, state, out));
}

Scalar compute_energy(Scalar TimeStep, const Energies& energies, const PhysicsState& state, const PhysicsState& state0) {
    // This function is responsible of computing the energy and derivatives for each energy in the simulation.
    // Here the energies must place the energy, gradient and Hessians to the correct place in the global energy and derivative structure
    Scalar phi = 0;
    energies.inertial_energies.for_each(ComputeInertialEnergyVisitor(TimeStep, state, state0, phi));
    energies.potential_energies.for_each(ComputePotentialEnergyVisitor(TimeStep, state, phi));
    return phi;
}

Vec compute_energy_gradient(Scalar TimeStep, const Energies& energies, const PhysicsState& state, const PhysicsState& state0) {
    // This function is responsible of computing the energy and derivatives for each energy in the simulation.
    // Here the energies must place the energy, gradient and Hessians to the correct place in the global energy and derivative structure
    Vec grad = Vec::Zero(state.get_nDoF());
    energies.inertial_energies.for_each(ComputeInertialEnergyGradientVisitor(TimeStep, state, state0, grad));
    energies.potential_energies.for_each(ComputePotentialEnergyGradientVisitor(TimeStep, state, grad));

    return grad;
}

void simulation_step(const Simulation& simulation, PhysicsState& state, EnergyAndDerivatives& f) {
    // Energy and derivatives computation
    const unsigned int nDoF = static_cast<unsigned int>(simulation.initial_state.x.size());

    // The state at the beginning of the step (x0, v0)
    const PhysicsState state0 = state;
    // Initial guess for our energy minimization
    // state =  PhysicsState(state0.x + simulation.TimeStep * state0.v, state0.v);

    const unsigned int maxIter = simulation.MaxNewtonIterations;
    // const Scalar reductionThreshold = 1e-4;
    // Scalar lastGradNorm = 1e9;
    for (unsigned int i = 0; i < maxIter; i++) {

        // Compute energy and derivatives
        // -----------------------------------------------------------------------------------------
        f.clear(nDoF);
        std::vector<ContactEvent> events;
        compute_energy_and_derivatives(simulation.TimeStep, simulation.energies, state, state0, f);
        // compute_energy_and_derivatives_finite(simulation.TimeStep, simulation.energies, state, state0, f);
        // const Scalar gradNorm = f.gradient.norm();
        // if (lastGradNorm - gradNorm < reductionThreshold) continue;
        // lastGradNorm = gradNorm;
        find_point_particle_contact_events(simulation.colliders, simulation.simulables, state, events);
        compute_contact_events_energy_and_derivatives(simulation.TimeStep, events, state, f);

        const Scalar energy0 = f.energy;

        // Integration step
        // -----------------------------------------------------------------------------------------
        Vec dx;
        integrate_implicit_euler(simulation, state, f, dx);

        // Update state
        // -----------------------------------------------------------------------------------------
        const PhysicsState stepState = state;
        update_simulation_state(simulation.TimeStep, simulation.energies, dx, state, state0);

        // Line search
        // -----------------------------------------------------------------------------------------
        const bool enableLineSearch = simulation.enable_line_search;
        if (not enableLineSearch) continue;
        Scalar lineSearchEnergy = compute_energy(simulation.TimeStep, simulation.energies, state, state0);
        DEBUG_LOG(energy0);
        DEBUG_LOG(lineSearchEnergy);
        Scalar alpha = 1.0f;
        const Scalar alpha_min_threshold = 1e-7;
        while (lineSearchEnergy > energy0) {
            DEBUG_LOG(alpha);
            if (alpha < alpha_min_threshold or std::isnan(energy0)) break;
            alpha /= 2.0f;
            state = stepState;
            update_simulation_state(simulation.TimeStep, simulation.energies, alpha * dx, state, state0);
            lineSearchEnergy = compute_energy(simulation.TimeStep, simulation.energies, state, state0);
        }
    }
}

void simulation_step(const Simulation& simulation, PhysicsState& state) {
    EnergyAndDerivatives f(0);
    simulation_step(simulation, state, f);
}

void update_simulation_state(const Scalar TimeStep, const Energies& energies, const Vec& dx, PhysicsState& state, const PhysicsState& state0) {
    energies.inertial_energies.for_each(UpdateSimulationStateVisitor(TimeStep, dx, state, state0));
}

inline Vec compute_energy_gradient_finite(Scalar E0, Scalar dx, Scalar TimeStep, const Energies& energies, const PhysicsState& state, const PhysicsState& state0) {
    const unsigned int nDoF = state.get_nDoF();
    Vec gradient = Vec::Zero(nDoF);
    for (unsigned int i = 0; i < nDoF; i++) {
        Vec dx_vec = Vec::Zero(nDoF);
        dx_vec[i] = dx;
        PhysicsState dstate = state;
        update_simulation_state(TimeStep, energies, dx_vec, dstate, state0);

        const Scalar dE = compute_energy(TimeStep, energies, dstate, state0);
        gradient[i] = (dE - E0) / dx;
    }
    return gradient;
}

inline Mat compute_energy_hessian_finite(Scalar E0, Scalar dx, Scalar TimeStep, const Energies& energies, const PhysicsState& state, const PhysicsState& state0) {
    const unsigned int nDoF = state.get_nDoF();
    Mat hessian = Mat::Zero(nDoF, nDoF);
    const Vec grad0 = compute_energy_gradient_finite(E0, dx, TimeStep, energies, state, state0);
    for (unsigned int i = 0; i < nDoF; i++) {
        Vec dx_vec = Vec::Zero(nDoF);
        dx_vec[i] = dx;
        PhysicsState dstate = state;
        update_simulation_state(TimeStep, energies, dx_vec, dstate, state0);

        const Scalar dE = compute_energy(TimeStep, energies, dstate, state0);
        const Vec grad = compute_energy_gradient_finite(dE, dx, TimeStep, energies, dstate, state0);
        hessian.row(i) = (grad - grad0) / dx;
    }
    return hessian;
}

void compute_energy_and_derivatives_finite(Scalar TimeStep, const Energies& energies, const PhysicsState& state, const PhysicsState& state0, EnergyAndDerivatives& out) {
    const Scalar dx = 1e-7;
    const unsigned int nDoF = state.get_nDoF();

    out.energy = compute_energy(TimeStep, energies, state, state0);
    out.gradient = compute_energy_gradient_finite(out.energy, dx, TimeStep, energies, state, state0);
    Mat hess = compute_energy_hessian_finite(out.energy, dx, TimeStep, energies, state, state0);
    const Scalar tol = dx;
    for (unsigned int i = 0; i < nDoF; i++) {
        for (unsigned int j = 0; j < nDoF; j++) {
            if (std::fabs(hess(i,j)) > tol) {
                out.hessian_triplets.emplace_back(i,j,hess(i,j));
            }
        }
    }
}
