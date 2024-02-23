#include <algorithm>
#include "linear_algebra.hpp"
#include "particle_rigid_body_copuling.hpp"
#include "utility_functions.hpp"

inline Eigen::Matrix<Scalar,3,6> compute_A_matrix(const ParticleRigidBodyCopuling& copuling, const PhysicsState& state) {
  const Mat3 skew_R0p = skew(copuling.rb.compute_rotation_matrix(state.x) * copuling.pos);
  Eigen::Matrix<Scalar,3,6> A = Eigen::Matrix<Scalar,3,6>::Zero();
  A.block<3,3>(0,0) = Mat3::Identity();
  A.block<3,3>(0,3) = - skew_R0p;
  return A;
}

void Copulings::compute_rigid_body_index_conversion() {
  for (unsigned int i = 0; i < copulings.size(); i++) {
    const ParticleRigidBodyCopuling& copuling = copulings[i];
    if (rigid_body_indices_conversion.contains(copuling.rb.index)) continue;
    rigid_body_indices_conversion[copuling.rb.index] = copuling.rb.index;
    for (unsigned int j = i; j < copulings.size(); j++) {
      const ParticleRigidBodyCopuling& other = copulings[j];
      // When a particle is avobe a rigid body, it will displace it's index by 3
      if (copuling.rb.index > other.particle.index) {
        rigid_body_indices_conversion[copuling.rb.index] -= 3;
      }
    }
  }
}

void Copulings::compute_dof_index_to_copuling() {
  for (unsigned int i = 0; i < copulings.size(); i++) {
    const ParticleRigidBodyCopuling& copuling = copulings[i];
    dof_index_to_copuling[copuling.particle.index] = i;
  }
}

void compute_copuling_jacobian(const Copulings& copulings, const PhysicsState& state, SparseMat& copuling_jacobian) {
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
    if (copulings.dof_index_to_copuling.contains(dof_index)) {
      const ParticleRigidBodyCopuling& copuling = copulings.get_copuling(dof_index);
      if (dof_index == copuling.particle.index) {
        index_offset+=3;
        const Eigen::Matrix<Scalar, 3, 6> A = compute_A_matrix(copuling, state);
        for (unsigned int j=0; j < 3; j++) {
          for (unsigned int k=0; k < 6; k++) {
            sparse_triplets.emplace_back(dof_index + j, copulings.rigid_body_indices_conversion.at(copuling.rb.index) + k, A(j,k));
          }
        }
        dof_index+=2; // skip the rest of particle's dof
      }
    }
    else {
      sparse_triplets.emplace_back(dof_index, dof_index - index_offset, 1);
    }
  }

  SparseMat fix_particle_jacobian(nDoF, nDoF-3*copulings.copulings.size()); // Reduces the dimensionality of the problem 3 dofs per copuled particle
  fix_particle_jacobian.setFromTriplets(sparse_triplets.begin(), sparse_triplets.end());
  copuling_jacobian = fix_particle_jacobian;
}
