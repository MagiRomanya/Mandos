#include <Eigen/Dense> // For inverse matrix

#include "fem_unit.hpp"
#include "utility_functions.hpp"

// F == deformation tensor
// epsilon == strain tensor
// sigma == stress tensor

template<unsigned int D>
Eigen::Vector<Scalar, D*D> vectorize_matrix(const Eigen::Matrix<Scalar,D,D>& mat) {
  Eigen::Vector<Scalar, D*D> result = Eigen::Vector<Scalar, D*D>::Zero();
  for (unsigned int i = 0; i< D; i++) {
    for (unsigned int j = 0; j< D; j++) {
      result(D*i+j) = mat(i,j);
    }
  }
  return result;
}

template<unsigned int N, unsigned int M>
Eigen::Matrix<Scalar, 3*N, 3*M> block_matrix(const Eigen::Matrix<Scalar,N,M>& mat) {
  Eigen::Matrix<Scalar, 3*N, 3*M> result = Eigen::Matrix<Scalar, 3*N, 3*M>::Zero();
  for (unsigned int i = 0; i< N; i++) {
    for (unsigned int j = 0; j< M; j++) {
      // result.block(3,3,i*3,j*3) = mat(i, j) * Mat3::Identity();
      result.template block<3,3>(i*3, j*3) = mat(i, j) * Mat3::Identity();
    }
  }
  return result;
}

Mat3 compute_strain_tensor_from_deformation_tensor(const Mat3& F) {
  return 0.5 * (F + F.transpose()) - Mat3::Identity();
}

Mat3 FEM_ElementParameters::compute_stress_tensor(const Mat3& epsilon) const {
  // return mu * ((F + F.transpose()) - 2 * I) + lambda *(F-I).trace() * I;
  return 2 * mu * epsilon + lambda * epsilon.trace() * Mat3::Identity();
}

Scalar FEM_ElementParameters::get_energy_density(const Mat3& epsilon) const {
  return mu * (epsilon.transpose()*epsilon).trace() + 0.5 * lambda * epsilon.trace() * epsilon.trace();
}

Eigen::Vector<Scalar, 12> FEM_ElementParameters::get_force_density(const Mat3& sigma) const {
  // NOTE We have to transpose the ds_dx matrix before making it a block matrix
  return - vectorize_matrix<3>(sigma).transpose() * block_matrix<3,4>(ds_dx.transpose());
}

Eigen::Matrix<Scalar, 12, 12> FEM_ElementParameters::get_df_dx_density() const {
  Eigen::Matrix<Scalar,9,9> d_vec_sigma_d_vec_F;
  d_vec_sigma_d_vec_F << 2*mu+lambda,  0,  0,  0,      lambda,  0,  0,  0,      lambda,
                         0          , mu,  0, mu,           0,  0,  0,  0,           0,
                         0          ,  0, mu,  0,           0,  0, mu,  0,           0,
                         0          , mu,  0, mu,           0,  0,  0,  0,           0,
                         lambda     ,  0,  0,  0, 2*mu+lambda,  0,  0,  0,      lambda,
                         0          ,  0,  0,  0,           0, mu,  0, mu,           0,
                         0          ,  0, mu,  0,           0,  0, mu,  0,           0,
                         0          ,  0,  0,  0,           0, mu,  0, mu,           0,
                         lambda     ,  0,  0,  0,      lambda,  0,  0,  0, 2*mu+lambda;

  const auto d_vec_F_d_x = block_matrix<3,4>(ds_dx.transpose());
  return d_vec_F_d_x.transpose() * d_vec_sigma_d_vec_F * d_vec_F_d_x;
}

Scalar FEM_ElementParameters::compute_volume(const Vec3& x1, const Vec3& x2, const Vec3& x3, const Vec3& x4) const {
  return compute_tetrahedron_volume(x2-x1, x3-x1, x4-x1);
}


Eigen::Matrix<Scalar,4,3> compute_shape_function_derivative(const Vec3& x1, const Vec3& x2, const Vec3& x3, const Vec3& x4) {
  Eigen::Matrix<Scalar,4,3> ds_dp;
  ds_dp << -1,-1,-1,
            1, 0, 0,
            0, 1, 0,
            0, 0, 1;

  Mat3 B;
  B << x2-x1, x3-x1, x4-x1;
  return ds_dp * B.inverse();
}
