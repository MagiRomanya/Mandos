#ifndef MANDOS_PARTICLE_RIGID_BODY_COUPLING_H_
#define MANDOS_PARTICLE_RIGID_BODY_COUPLING_H_

#include <vector>

#include <Mandos/linear_algebra.hpp>
#include <Mandos/particle.hpp>
#include <Mandos/physics_state.hpp>
#include <Mandos/rigid_body.hpp>

namespace mandos
{
struct ParticleRigidBodyCoupling {
    ParticleRigidBodyCoupling(const RigidBody& rb, const Particle& p, Vec3 rel_pos)
        : rb(rb)
        , particle(p)
        , pos(rel_pos)
    {
    }

    const RigidBody rb;
    const Particle particle;
    const Vec3 pos;  // Particle fixed position wrt the RB rotating frame
};

struct Couplings {
    Couplings(){};
    Couplings(const std::vector<ParticleRigidBodyCoupling>& couplings)
        : couplings(couplings)
    {
        compute_rigid_body_index_conversion();
        compute_dof_index_to_coupling();
    }

    inline void add_coupling(ParticleRigidBodyCoupling coupling)
    {
        couplings.push_back(coupling);
        compute_rigid_body_index_conversion();
        compute_dof_index_to_coupling();
    }

    std::vector<ParticleRigidBodyCoupling> couplings;

    // Rigid body index (full dof state) --> Rigid body index (copuled dofs state)
    std::unordered_map<unsigned int, unsigned int> rigid_body_indices_conversion;

    // Dof index (full dof state) --> coupling index
    std::unordered_map<unsigned int, unsigned int> dof_index_to_coupling;

    void compute_rigid_body_index_conversion();
    void compute_dof_index_to_coupling();

    inline ParticleRigidBodyCoupling get_coupling(unsigned int dof_index) const
    {
        return couplings[dof_index_to_coupling.at(dof_index)];
    }
};

void compute_coupling_jacobian(const Couplings& couplings, const PhysicsState& state, SparseMat& coupling_jacobian);

}  // namespace mandos

#endif  // MANDOS_PARTICLE_RIGID_BODY_COPULING_H_
