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


struct Copulings {
    Copulings() {};
    Copulings(const std::vector<ParticleRigidBodyCopuling>& copulings) : copulings(copulings) {
        compute_rigid_body_index_conversion();
        compute_dof_index_to_copuling();
    }

    std::vector<ParticleRigidBodyCopuling> copulings;

    // Rigid body index (full dof state) --> Rigid body index (copuled dofs state)
    std::unordered_map<unsigned int, unsigned int> rigid_body_indices_conversion;

    // Dof index (full dof state) --> copuling index
    std::unordered_map<unsigned int, unsigned int> dof_index_to_copuling;

    void compute_rigid_body_index_conversion();
    void compute_dof_index_to_copuling();

    inline ParticleRigidBodyCopuling get_copuling(unsigned int dof_index) const {return copulings[dof_index_to_copuling.at(dof_index)];}

};

void compute_copuling_jacobian(const Copulings& copulings, const PhysicsState& state, SparseMat& copuling_jacobian);


#endif // PARTICLE_RIGID_BODY_COPULING_H_
