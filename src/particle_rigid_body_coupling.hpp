#ifndef PARTICLE_RIGID_BODY_COUPLING_H_
#define PARTICLE_RIGID_BODY_COUPLING_H_

#include "linear_algebra.hpp"
#include "particle.hpp"
#include "physics_state.hpp"
#include "rigid_body.hpp"
#include <vector>

struct ParticleRigidBodyCouplings {
    ParticleRigidBodyCouplings(const RigidBody& rb, const Particle& p, Vec3 rel_pos)
        : rb(rb), particle(p), pos(rel_pos) {}

    const RigidBody rb;
    const Particle particle;
    const Vec3 pos;                         // Particle fixed position wrt the RB rotating frame
};


struct Couplings {
    Couplings() {};
    Couplings(const std::vector<ParticleRigidBodyCouplings>& couplings) : couplings(couplings) {
        compute_rigid_body_index_conversion();
        compute_dof_index_to_coupling();
    }

    inline void add_coupling(ParticleRigidBodyCouplings coupling) {
        couplings.push_back(coupling);
        compute_rigid_body_index_conversion();
        compute_dof_index_to_coupling();
    }

    std::vector<ParticleRigidBodyCouplings> couplings;

    // Rigid body index (full dof state) --> Rigid body index (coupled dofs state)
    std::unordered_map<unsigned int, unsigned int> rigid_body_indices_conversion;

    // Dof index (full dof state) --> coupling index
    std::unordered_map<unsigned int, unsigned int> dof_index_to_coupling;

    void compute_rigid_body_index_conversion();
    void compute_dof_index_to_coupling();

    inline ParticleRigidBodyCouplings get_couplings(unsigned int dof_index) const {return couplings[dof_index_to_coupling.at(dof_index)];}
};

void compute_coupling_jacobian(const Couplings& couplings, const PhysicsState& state, SparseMat& coupling_jacobian);


#endif // PARTICLE_RIGID_BODY_COUPLING_H_
