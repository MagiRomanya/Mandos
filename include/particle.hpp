#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "linear_algebra.hpp"
#include "physics_state.hpp"

struct Particle {
    Particle(Scalar mass, unsigned int index) : index(index), mass(mass) {}
    const unsigned int index;
    const Scalar mass;

    Vec3 get_position(const PhysicsState& state) const { return state.x.segment(index, 3); }
    Vec3 get_velocity(const PhysicsState& state) const { return state.v.segment(index, 3); }

    // Scalar get_kinetic_energy(const PhysicsState& state) const { return 0.5 * mass * get_velocity(state).squaredNorm(); }
};

#endif // PARTICLE_H_
