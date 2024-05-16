#ifndef MANDOS_PARTICLE_H_
#define MANDOS_PARTICLE_H_

#include <Mandos/linear_algebra.hpp>
#include <Mandos/physics_state.hpp>

namespace mandos
{

struct Particle {
    Particle(Scalar mass, unsigned int index)
        : index(index)
        , mass(mass)
    {
    }

    unsigned int index;
    Scalar mass;

    inline Vec3 get_position(const PhysicsState& state) const
    {
        return state.x.segment<3>(index);
    }

    inline Vec3 get_velocity(const PhysicsState& state) const
    {
        return state.v.segment<3>(index);
    }
};

}  // namespace mandos

#endif  // MANDOS_PARTICLE_H_
