#include <algorithm>

#include <Mandos/linear_algebra.hpp>
#include <Mandos/particle_rigid_body_coupling.hpp>
#include <Mandos/utility_functions.hpp>

namespace mandos
{

inline Eigen::Matrix<Scalar, 3, 6> compute_A_matrix(const ParticleRigidBodyCoupling& coupling,
                                                    const PhysicsState& state)
{
    const Mat3 skew_R0p = skew(coupling.rb.compute_rotation_matrix(state.x) * coupling.pos);
    Eigen::Matrix<Scalar, 3, 6> A = Eigen::Matrix<Scalar, 3, 6>::Zero();
    A.block<3, 3>(0, 0) = Mat3::Identity();
    A.block<3, 3>(0, 3) = -skew_R0p;
    return A;
}

void Couplings::compute_rigid_body_index_conversion()
{
    for (unsigned int i = 0; i < couplings.size(); i++) {
        const ParticleRigidBodyCoupling& coupling = couplings[i];
        if (rigid_body_indices_conversion.contains(coupling.rb.index))
            continue;
        rigid_body_indices_conversion[coupling.rb.index] = coupling.rb.index;
        for (unsigned int j = i; j < couplings.size(); j++) {
            const ParticleRigidBodyCoupling& other = couplings[j];
            // When a particle is avobe a rigid body, it will displace it's index by 3
            if (coupling.rb.index > other.particle.index) {
                rigid_body_indices_conversion[coupling.rb.index] -= 3;
            }
        }
    }
}

void Couplings::compute_dof_index_to_coupling()
{
    for (unsigned int i = 0; i < couplings.size(); i++) {
        const ParticleRigidBodyCoupling& coupling = couplings[i];
        dof_index_to_coupling[coupling.particle.index] = i;
    }
}

void compute_coupling_jacobian(const Couplings& couplings, const PhysicsState& state, SparseMat& coupling_jacobian)
{
    const unsigned int nDoF = state.get_nDoF();

    std::vector<Triplet> sparse_triplets;

    /**
     * Construct the fixed particle jacobian matrix
     * dp = A rb
     * / dp  \   /A  0\ / dp  \
     * | rb  | = |I  0| | rb  |
     * \other/   \0  I/ \other/
     */
    unsigned int index_offset = 0;
    for (unsigned int dof_index = 0; dof_index < nDoF; dof_index++) {
        if (couplings.dof_index_to_coupling.contains(dof_index)) {
            const ParticleRigidBodyCoupling& coupling = couplings.get_coupling(dof_index);
            if (dof_index == coupling.particle.index) {
                index_offset += 3;
                const Eigen::Matrix<Scalar, 3, 6> A = compute_A_matrix(coupling, state);
                for (unsigned int j = 0; j < 3; j++) {
                    for (unsigned int k = 0; k < 6; k++) {
                        sparse_triplets.emplace_back(
                            dof_index + j, couplings.rigid_body_indices_conversion.at(coupling.rb.index) + k, A(j, k));
                    }
                }
                dof_index += 2;  // skip the rest of particle's dof
            }
        } else {
            sparse_triplets.emplace_back(dof_index, dof_index - index_offset, 1);
        }
    }

    SparseMat fix_particle_jacobian(
        nDoF,
        nDoF -
            3 * couplings.couplings.size());  // Reduces the dimensionality of the problem 3 dofs per copuled particle
    fix_particle_jacobian.setFromTriplets(sparse_triplets.begin(), sparse_triplets.end());
    coupling_jacobian = fix_particle_jacobian;
}

}  // namespace mandos