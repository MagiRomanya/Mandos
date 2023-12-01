#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "linear_algebra.hpp"
#include "physics_state.hpp"
#include <iostream>

struct Particle {
    Particle(Scalar mass, unsigned int index) : index(index), mass(mass) {}
    const unsigned int index;
    const Scalar mass;

    inline Vec3 get_position(const PhysicsState& state) const { return state.x.segment(index, 3); }
    inline Vec3 get_position_old(const PhysicsState& state) const { return state.x_old.segment(index, 3); }
    inline Vec3 get_position_old2(const PhysicsState& state) const { return state.x_old2.segment(index, 3); }

    inline Vec3 get_velocity(const PhysicsState& state, Scalar TimeStep) const {
        return (get_position(state) - get_position_old(state)) / TimeStep;
    }
    inline Vec3 get_velocity_old(const PhysicsState& state, Scalar TimeStep) const {
        return (get_position_old(state) - get_position_old2(state)) / TimeStep;
    }

    inline void update_state(const Vec& dx, Vec& x) const {
        x.segment<3>(index) += dx.segment<3>(index);
    }
};

#endif // PARTICLE_H_
