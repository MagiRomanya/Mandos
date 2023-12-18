#ifndef PARTICLE_RIGID_BODY_COPULING_H_
#define PARTICLE_RIGID_BODY_COPULING_H_

#include "linear_algebra.hpp"
#include "particle.hpp"
#include "physics_state.hpp"
#include "rigid_body.hpp"
#include <vector>

struct ParticleRigidBodyCopuling {
    ParticleRigidBodyCopuling(const RigidBody& rb, const Particle& p, Vec3 rel_pos)
        : rb(rb), particle(p), pos(rel_pos) {}

    const RigidBody rb;
    const Particle particle;
    const Vec3 pos;                         // Particle fixed position wrt the RB rotating frame

};

bool operator<(const ParticleRigidBodyCopuling& c1, const ParticleRigidBodyCopuling& c2);

void compute_copuling_jacobian(const std::vector<ParticleRigidBodyCopuling> copulings, const PhysicsState& state, SparseMat& copuling_jacobian);


#endif // PARTICLE_RIGID_BODY_COPULING_H_
