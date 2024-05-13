#pragma once

#include "fem_element.hpp"
#include "gravity.hpp"
#include "inertia_energies.hpp"
#include "spring.hpp"
#include "rod_segment.hpp"
#include "utility_functions.hpp"

struct InertialEnergies {
    std::vector<LinearInertia> linear_inertias;
    std::vector<RotationalInertia> rotational_inertias;
    std::vector<RotationalInertiaGlobal> rotational_global_inertias;

    template<typename Visitor>
    void for_each(Visitor visitor) const {
        visitor(linear_inertias);
        visitor(rotational_inertias);
        visitor(rotational_global_inertias);
    }
};

struct PotentialEnergies {
    std ::vector<ParticleSpring> particle_springs;
    std ::vector<RigidBodySpring> rigid_body_springs;
    std ::vector<RodSegment> rod_segments;
    std ::vector<FEM_Element3D<FEM_LinearMaterial>> fem_elements_linearMat;
    std ::vector<FEM_Element3D<FEM_NeoHookeanMaterial>> fem_elements_neoHookMat;
    std ::vector<Gravity> gravities;

    template<typename Visitor>
    void for_each(Visitor visitor) const {
        visitor(gravities);
        visitor(particle_springs);
        visitor(rigid_body_springs);
        visitor(rod_segments);
        visitor(fem_elements_linearMat);
        visitor(fem_elements_neoHookMat);
    }
};

struct Energies {
    InertialEnergies inertial_energies;
    PotentialEnergies potential_energies;
};

template <typename MaterialType>
void add_FEM_element(Energies& energies, FEM_Element3D<MaterialType> element);

template void add_FEM_element(Energies &energies, FEM_Element3D<FEM_LinearMaterial> element);
template void add_FEM_element(Energies &energies, FEM_Element3D<FEM_NeoHookeanMaterial> element);

template <typename MaterialType>
inline void add_FEM_element(Energies& energies, FEM_Element3D<MaterialType> element) {
    if constexpr (std ::is_same<MaterialType, FEM_LinearMaterial>()) {
        energies.potential_energies.fem_elements_linearMat.push_back(element);
    }
    if constexpr (std ::is_same<MaterialType, FEM_NeoHookeanMaterial>()) {
        energies.potential_energies.fem_elements_neoHookMat.push_back(element);
    }
}

//////////////////////////////////////////////////////////////////////////////
////////////////////////////// ENERGY VISITORS////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

struct ComputeInertialEnergyVisitor {
    ComputeInertialEnergyVisitor(Scalar TimeStep, const PhysicsState& state, const PhysicsState& state0, Scalar& out)
        : TimeStep(TimeStep), state(state), state0(state0), out(out) {}
    Scalar TimeStep;
    const PhysicsState& state;
    const PhysicsState& state0;
    Scalar& out;

    template <typename T>
    void operator()(const std::vector<T>& energies) {
        for (unsigned int i = 0; i < energies.size(); i++) {
            out += energies[i].compute_energy(TimeStep, state, state0);
        }
    }
};

struct ComputePotentialEnergyVisitor {
    ComputePotentialEnergyVisitor(Scalar TimeStep, const PhysicsState& state, Scalar& out)
        : TimeStep(TimeStep), state(state), out(out) {}
    Scalar TimeStep;
    const PhysicsState& state;
    Scalar& out;

    template <typename T>
    void operator()(const std::vector<T>& energies) {
        for (unsigned int i = 0; i < energies.size(); i++) {
            out += energies[i].compute_energy(TimeStep, state);
        }
    }
};

struct ComputeInertialEnergyGradientVisitor {
    ComputeInertialEnergyGradientVisitor(Scalar TimeStep, const PhysicsState& state, const PhysicsState& state0, Vec& out)
        : TimeStep(TimeStep), state(state), state0(state0), out(out) {}
    Scalar TimeStep;
    const PhysicsState& state;
    const PhysicsState& state0;
    Vec& out;

    template <typename T>
    void operator()(const std::vector<T>& energies) {
        for (unsigned int i = 0; i < energies.size(); i++) {
            energies[i].compute_energy_gradient(TimeStep, state, state0, out);
        }
    }
};

struct ComputePotentialEnergyGradientVisitor {
    ComputePotentialEnergyGradientVisitor(Scalar TimeStep, const PhysicsState& state, Vec& out)
        : TimeStep(TimeStep), state(state), out(out) {}
    Scalar TimeStep;
    const PhysicsState& state;
    Vec& out;

    template <typename T>
    void operator()(const std::vector<T>& energies) {
        for (unsigned int i = 0; i < energies.size(); i++) {
            energies[i].compute_energy_gradient(TimeStep, state, out);
        }
    }
};

struct ComputeInertialEnergyAndDerivativesVisitor {
    ComputeInertialEnergyAndDerivativesVisitor(Scalar TimeStep, const PhysicsState& state, const PhysicsState& state0, EnergyAndDerivatives& out)
        : TimeStep(TimeStep), state(state), state0(state0), out(out) {}
    Scalar TimeStep;
    const PhysicsState& state;
    const PhysicsState& state0;
    EnergyAndDerivatives& out;

    template <typename T>
    void operator()(const std::vector<T>& energies) {
        for (unsigned int i = 0; i < energies.size(); i++) {
            energies[i].compute_energy_and_derivatives(TimeStep, state, state0, out);
        }
    }
};

struct ComputePotentialEnergyAndDerivativesVisitor {
    ComputePotentialEnergyAndDerivativesVisitor(Scalar TimeStep, const PhysicsState& state, EnergyAndDerivatives& out)
        : TimeStep(TimeStep), state(state), out(out) {}
    Scalar TimeStep;
    const PhysicsState& state;
    EnergyAndDerivatives& out;

    template <typename T>
    void operator()(const std::vector<T>& energies) {
        for (unsigned int i = 0; i < energies.size(); i++) {
            energies[i].compute_energy_and_derivatives(TimeStep, state, out);
            // DEBUG_LOG(out.energy);
        }
    }
};

struct UpdateSimulationStateVisitor {
    UpdateSimulationStateVisitor(Scalar TimeStep, const Vec& dx, PhysicsState& state, const PhysicsState& state0)
        : TimeStep(TimeStep), dx(dx), state(state), state0(state0) {}
    Scalar TimeStep;
    const Vec& dx;
    PhysicsState& state;
    const PhysicsState& state0;

    template <typename T>
    void operator()(const std::vector<T>& energies) {
        for (unsigned int i = 0; i < energies.size(); i++) {
            energies[i].update_state(TimeStep, dx, state, state0);
        }
    }
};
