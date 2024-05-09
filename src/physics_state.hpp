#ifndef PHYSICS_STATE_H_
#define PHYSICS_STATE_H_

#include <vector>
#include "linear_algebra.hpp"

struct PhysicsState {
    PhysicsState() {}
    PhysicsState(const Vec& x, const Vec& v) : x(x), v(v) {}
    Vec x;
    Vec v;

    inline void add_size(unsigned int increment_dof) {
        x.conservativeResize(x.size() + increment_dof);
        v.conservativeResize(v.size() + increment_dof);
    }

    inline unsigned int get_nDoF() const {return static_cast<unsigned int>( x.size() ); }
};

struct EnergyAndDerivatives {
    EnergyAndDerivatives() {}
    EnergyAndDerivatives(unsigned int nDoF) {
        energy = 0;
        gradient.setZero(nDoF);
        hessian_triplets.clear();
    }
    inline void clear(unsigned int nDoF) {
        energy = 0;
        gradient.setZero(nDoF);
        hessian_triplets.clear();
    }
    // Container
    Scalar energy;
    Vec gradient;
    std::vector<Triplet> hessian_triplets;
};

/**
 * Inertial Energy
 *
 * Virtual class exposing the mandatory API that each inertia energy must implement.
 * Note that this virtual class is never instantiated at runtime, and only exists
 * for nicer compilation time errors.
 */
struct InertialEnergy {
    virtual Scalar compute_energy(Scalar TimeStep, const PhysicsState& state, const PhysicsState& state0) const = 0;
    virtual void compute_energy_gradient(Scalar TimeStep, const PhysicsState& state, const PhysicsState& state0, Vec& grad) const = 0;
    virtual void compute_energy_and_derivatives(Scalar TimeStep, const PhysicsState& state, const PhysicsState& state0, EnergyAndDerivatives& out) const = 0;

    virtual void update_state(const Scalar TimeStep, const Vec& dx, PhysicsState& state, const PhysicsState& state0) const = 0;
};

/**
 * Potential Energy
 *
 * Virtual class exposing the mandatory API that each potential energy must implement.
 * Note that this virtual class is never instantiated at runtime, and only exists
 * for nicer compilation time errors.
 */
struct PotentialEnergy {
    virtual Scalar compute_energy(Scalar TimeStep, const PhysicsState& state) const = 0;
    virtual void compute_energy_gradient(Scalar TimeStep, const PhysicsState& state, Vec& grad) const = 0;
    virtual void compute_energy_and_derivatives(Scalar TimeStep, const PhysicsState& state, EnergyAndDerivatives& out) const = 0;
};

#endif // PHYSICS_STATE_H_
