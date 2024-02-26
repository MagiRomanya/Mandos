#include <cmath>
#include <vector>

#include "fem_unit.hpp"
#include "simulation.hpp"
#include "integrators.hpp"
#include "linear_algebra.hpp"
#include "particle.hpp"
#include "physics_state.hpp"
#include "rigid_body.hpp"
#include "utility_functions.hpp"

void compute_energy_and_derivatives_finite(Scalar TimeStep, const Energies& energies, const PhysicsState& state, const PhysicsState& state0, EnergyAndDerivatives& out);

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

    // The state at the beginning of the step (x0, v0)
    const PhysicsState state0 = state;
    // Initial guess for our energy minimization
    // state =  PhysicsState(state0.x + simulation.TimeStep * state0.v, state0.v);

    const unsigned int maxIter = 3;
    for (unsigned int i = 0; i < maxIter; i++) {

        // Compute energy and derivatives
        // -----------------------------------------------------------------------------------------
        f = EnergyAndDerivatives(nDoF);
        compute_energy_and_derivatives(simulation.TimeStep, simulation.energies, state, state0, f);
        const Scalar energy0 = f.energy;
        // DEBUG_LOG(f.gradient.norm());

        // Integration step
        // -----------------------------------------------------------------------------------------
        Vec dx;
        integrate_implicit_euler(simulation, state, f, dx);

        // Update state
        // -----------------------------------------------------------------------------------------
        const PhysicsState stepState = state;
        update_simulation_state(simulation.simulables, dx, state.x); // x_new
        state.v = (state.x - state0.x) / simulation.TimeStep; // v_new

        // Line search
        // -----------------------------------------------------------------------------------------
        const bool enableLineSearch = false;
        if (not enableLineSearch) continue;
        Scalar lineSearchEnergy = compute_energy(simulation.TimeStep, simulation.energies, state, state0);
        Scalar alpha = 1.0f;
        const float alpha_min_threshold = 1e-7;
        // DEBUG_LOG(energy);
        // DEBUG_LOG(lineSearchEnergy);
        while (lineSearchEnergy > energy0) {
            if (alpha < alpha_min_threshold or std::isnan(energy0)) break;
            alpha /= 2.0f;
            state = stepState;
            update_simulation_state(simulation.simulables, alpha * dx, state.x);
            state.v = (state.x - state0.x) / simulation.TimeStep; // v_new
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

inline Vec compute_energy_gradient_finite(Scalar E0, Scalar dx, Scalar TimeStep, const Energies& energies, const PhysicsState& state, const PhysicsState& state0) {
    const unsigned int nDoF = state.get_nDoF();
    Vec gradient = Vec::Zero(nDoF);
    for (unsigned int i = 0; i < nDoF; i++) {
        PhysicsState dstate = state;
        dstate.x[i] += dx;
        dstate.v[i] += dx / TimeStep;

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
        PhysicsState dstate = state;
        dstate.x[i] += dx;
        dstate.v[i] += dx / TimeStep;

        const Scalar dE = compute_energy(TimeStep, energies, dstate, state0);
        const Vec grad = compute_energy_gradient_finite(dE, dx, TimeStep, energies, dstate, state0);
        hessian.row(i) = (grad - grad0) / dx;
    }
    return hessian;
}

void compute_energy_and_derivatives_finite(Scalar TimeStep, const Energies& energies, const PhysicsState& state, const PhysicsState& state0, EnergyAndDerivatives& out) {
    const Scalar dx = 1e-8;
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
