#include <algorithm>
#include "linear_algebra.hpp"
#include "particle_rigid_body_copuling.hpp"
#include "utility_functions.hpp"

inline void compute_row_bounds(std::vector<Triplet>& row_wise_triplets, unsigned int index, unsigned int& start, unsigned int& end) {
  unsigned int width = 0;
  int last_row = -10;
  for (unsigned int i = 0; i< row_wise_triplets.size(); i++) {
    const Triplet& t = row_wise_triplets[i];

    if (last_row < (int)index) {
      if (t.row() == (int)index) start = i;
      else if (t.row() > (int)(index + width)) {
        start = i;
        end = i;
        break;
      }
    }
    else if ((last_row == (int)(index + width)) && (t.row() > (int)(index + width))) {
      assert(i==0);
      end = i-1;
      break;
    }
    last_row = t.row();
  }
  if (end == 0 && start !=0) end = row_wise_triplets.size()-1;
}

void ParticleRigidBodyCopuling::apply_copuling(const PhysicsState& state, SparseMat& hessian, Vec& gradient) const {
  const unsigned int nDoF = state.get_nDoF();
  const Mat3 skew_R0p = skew(rb.compute_rotation_matrix(state.x) * pos);
  Eigen::Matrix<Scalar,3,6> A = Eigen::Matrix<Scalar,3,6>::Zero();
  A.block<3,3>(0,0) = Mat3::Identity();
  A.block<3,3>(0,3) = - skew_R0p;
  std::vector<Triplet> sparse_triplets;

  // Construct the fixed particle jacobian matrix
  // dp = A rb
  // / dp  \   /A  0\ / rb  \
  // | rb  | = |I  0| |     |
  // \other/   \0  I/ \other/
  unsigned int index_offset = 0;
  for (unsigned int i = 0; i < nDoF; i++) {
    if (i == particle.index) {
      index_offset+=3;
      for (unsigned int j=0; j < 3; j++) {
        for (unsigned int k=0; k < 6; k++) {
          sparse_triplets.emplace_back(i, rb.index - index_offset, A(j,k));
        }
      }
      i+=3;
    }
    else {
      sparse_triplets.emplace_back(i, i - index_offset, 1);
    }
  }

  SparseMat fix_particle_jacobian(nDoF, nDoF-3); // Reduces the dimensionality of the problem by one particle
  fix_particle_jacobian.setFromTriplets(sparse_triplets.begin(), sparse_triplets.end());

  // Apply the jacobian to the global hessian and gradient
  const SparseMat fix_particle_jacobian_t = fix_particle_jacobian.transpose();
  hessian = fix_particle_jacobian_t * hessian * fix_particle_jacobian;
  gradient = fix_particle_jacobian_t * gradient;
}

// void ParticleRigidBodyCopuling::apply_copuling(const PhysicsState& state, std::vector<Triplet>& sorted_hessian_triplets, Vec& gradient) const {
//   gradient.segment<3>(particle.index) = Vec3::Zero();

//   const Mat3 skew_R0p = skew(rb.compute_rotation_matrix(state.x) * pos);

//   std::vector<Triplet> hess;
//   hess.reserve(sorted_hessian_triplets.size());
//   unsigned int particle_hessian_start_index = 0;
//   unsigned int particle_hessian_end_index = 0;
//   compute_row_bounds(sorted_hessian_triplets, particle.index, particle_hessian_start_index, particle_hessian_end_index);

//   for (unsigned int i = 0; i < particle_hessian_start_index; i++)
//     hess.push_back(sorted_hessian_triplets[i]);

//   // Construct the hessian of the particle respecting triplet row-wise ordering
//   for (unsigned int i = 0; i < 3; i++) {
//     if (particle.index < rb.index)
//       hess.emplace_back(particle.index+i, particle.index+i, 1);           // I

//     hess.emplace_back(particle.index+i, rb.index+i, -1);                  // -I
//     hess.emplace_back(particle.index+i, rb.index+3+0, skew_R0p(0,0));     // -skew(R0 p)
//     hess.emplace_back(particle.index+i, rb.index+3+1, skew_R0p(0,0));
//     hess.emplace_back(particle.index+i, rb.index+3+2, skew_R0p(0,0));

//     if (particle.index > rb.index)
//       hess.emplace_back(particle.index+i, particle.index+i, 1);           // I
//   }
//   compute_row_bounds(sorted_hessian_triplets, particle.index+3, particle_hessian_start_index, particle_hessian_end_index);
//   for (unsigned int i = 0; i < particle_hessian_end_index; i++)
//     hess.push_back(sorted_hessian_triplets[i]);
// }

void sort_triplets_row_major(std::vector<Triplet>& triplets) {
  std::sort(triplets.begin(), triplets.end(),
            [](const Triplet& a, const Triplet& b ) {
              return ((a.row() != b.row()) ? (a.row() < b.row()) : (a.col() < b.col())); // Row-Major
            });
}

void sort_triplets_column_major(std::vector<Triplet>& triplets) {
  std::sort(triplets.begin(), triplets.end(),
            [](const Triplet& a, const Triplet& b ) {
              return ((a.col() != b.col()) ? (a.col() < b.col()) : (a.row() < b.row())); // Column-Major
            });
}
