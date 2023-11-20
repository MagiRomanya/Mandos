#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "linear_algebra.hpp"
#include "physics_state.hpp"

struct Particle {
    Particle(Scalar mass, unsigned int index) : index(index), mass(mass) {}
    const unsigned int index;
    const Scalar mass;

    inline Vec3 get_position(const PhysicsState& state) const { return state.x.segment(index, 3); }
    inline Vec3 get_velocity(const PhysicsState& state) const { return state.v.segment(index, 3); }
    inline Scalar get_kinetic_energy(const PhysicsState& state) const { return 0.5 * mass * get_velocity(state).squaredNorm(); }

    inline void compute_energy_and_derivatives(const PhysicsState& state, EnergyAndDerivatives& out) const {
        out.energy += get_kinetic_energy(state);
    }

    inline void update_state(Scalar TimeStep, const Vec& new_velocities, PhysicsState& state) const {
        state.v.segment(index, 3) = new_velocities.segment(index, 3);
        state.x.segment(index, 3) += TimeStep * new_velocities.segment(index, 3);
    }
};

#endif // PARTICLE_H_
