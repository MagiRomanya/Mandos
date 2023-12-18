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

bool operator<(const ParticleRigidBodyCopuling& c1, const ParticleRigidBodyCopuling& c2) {
  return c1.particle.index < c2.particle.index;
}

void compute_copuling_jacobian(const std::vector<ParticleRigidBodyCopuling> copulings, const PhysicsState& state, SparseMat& copuling_jacobian) {
  const unsigned int nDoF = state.get_nDoF();

  std::vector<Triplet> sparse_triplets;

  // Construct the fixed particle jacobian matrix
  // dp = A rb
  // / dp  \   /A  0\ / dp  \
  // | rb  | = |I  0| | rb  |
  // \other/   \0  I/ \other/
  unsigned int index_offset = 0;
  unsigned int copuling_index = 0;
  bool are_copulings_left = not copulings.empty();

  for (unsigned int i = 0; i < nDoF; i++) {
    if (are_copulings_left and (i == copulings[copuling_index].particle.index)) {
      const ParticleRigidBodyCopuling& copuling = copulings[copuling_index];
      if (i == copuling.particle.index) {
        index_offset+=3;
        const Eigen::Matrix<Scalar, 3, 6> A = compute_A_matrix(copuling, state);
        for (unsigned int j=0; j < 3; j++) {
          for (unsigned int k=0; k < 6; k++) {
            sparse_triplets.emplace_back(i + j, copuling.rb.index - index_offset + k, A(j,k));
          }
        }
        i+=2; // skip the rest of particle's dof
        copuling_index++; // Check the next copuling (we are assuming ordered copulings)
        are_copulings_left = copuling_index + 1 == copulings.size();
      }
    }
    else {
      sparse_triplets.emplace_back(i, i - index_offset, 1);
    }
  }

  SparseMat fix_particle_jacobian(nDoF, nDoF-3*copulings.size()); // Reduces the dimensionality of the problem by one particle
  fix_particle_jacobian.setFromTriplets(sparse_triplets.begin(), sparse_triplets.end());
  copuling_jacobian = fix_particle_jacobian;
}
