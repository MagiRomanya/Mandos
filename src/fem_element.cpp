#include <Eigen/Dense> // For inverse matrix
#include "fem_unit.hpp"
#include "utility_functions.hpp"

// F == deformation tensor
// epsilon == strain tensor
// sigma == stress tensor

// Given a square DxD matrix mat, return a vector of D² components
// A  B
// C  D  ->  vec = A, B, C, D
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

// Given a NxM matrix, return a 3Nx3M matrix where each matrix component m_i_j becomes now m_i_j*I
// where I is the 3x3 identity matrix
// A  B            A* I, B*I
// C  D  ->  mat = C* I, D*I
template<unsigned int N, unsigned int M>
Eigen::Matrix<Scalar, 3*N, 3*M> block_matrix(const Eigen::Matrix<Scalar,N,M>& mat) {
  Eigen::Matrix<Scalar, 3*N, 3*M> result = Eigen::Matrix<Scalar, 3*N, 3*M>::Zero();
  for (unsigned int i = 0; i< N; i++) {
    for (unsigned int j = 0; j< M; j++) {
      // Wierd syntax we have to use because of templates
      // It just means result.block...
      result.template block<3,3>(i*3, j*3) = mat(i, j) * Mat3::Identity();
    }
  }
  return result;
}


// NOTE We have to transpose the ds_dx matrix before making it a block matrix
FEM_ElementParameters::FEM_ElementParameters(Scalar mu, Scalar lambda, Eigen::Matrix<Scalar,4,3> ds_dx)
  : mu(mu), lambda(lambda), dvecF_dx(block_matrix<3,4>(ds_dx.transpose())) {}

Mat3 FEM_ElementParameters::compute_deformation_tensor(const Vec3& x1, const Vec3& x2, const Vec3& x3, const Vec3& x4) const {
  // X_vec is a 12 component column vector organized as X_vec = x1, x2, x3, x4
  Eigen::Vector<Scalar, 12> X_vec = Eigen::Vector<Scalar, 12>::Zero();
  X_vec.segment<3>(3*0) = x1;
  X_vec.segment<3>(3*1) = x2;
  X_vec.segment<3>(3*2) = x3;
  X_vec.segment<3>(3*3) = x4;
  // F_vec is a column-wise vectorized representation of the F tensor
  // Meaning F_vec = F11, F21, F31, F12, F22, F32, F13, F32, F33
  const Eigen::Vector<Scalar, 9> F_vec = dvecF_dx * X_vec;
  Mat3 F;
  F << F_vec(0), F_vec(3), F_vec(6),
       F_vec(1), F_vec(4), F_vec(7),
       F_vec(2), F_vec(5), F_vec(8);
  return F;
}

Mat3 compute_strain_tensor_from_deformation_tensor(const Mat3& F) {
  return 0.5 * (F + F.transpose()) - Mat3::Identity();
}

Mat3 FEM_ElementParameters::compute_stress_tensor(const Mat3& epsilon) const {
  // Homogenious material Hook's Law https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
  return 2 * mu * epsilon + lambda * epsilon.trace() * Mat3::Identity();
}

Scalar FEM_ElementParameters::get_energy_density(const Mat3& epsilon) const {
  return mu * (epsilon.transpose()*epsilon).trace() + 0.5 * lambda * epsilon.trace() * epsilon.trace();
}

Eigen::Vector<Scalar, 12> FEM_ElementParameters::get_force_density(const Mat3& sigma) const {
  return - vectorize_matrix<3>(sigma).transpose() * dvecF_dx;
}

Eigen::Matrix<Scalar, 12, 12> FEM_ElementParameters::get_df_dx_density() const {
  Eigen::Matrix<Scalar,9,9> dvecSigma_dvecF;
  dvecSigma_dvecF << 2*mu+lambda,  0,  0,  0,      lambda,  0,  0,  0,      lambda,
                     0          , mu,  0, mu,           0,  0,  0,  0,           0,
                     0          ,  0, mu,  0,           0,  0, mu,  0,           0,
                     0          , mu,  0, mu,           0,  0,  0,  0,           0,
                     lambda     ,  0,  0,  0, 2*mu+lambda,  0,  0,  0,      lambda,
                     0          ,  0,  0,  0,           0, mu,  0, mu,           0,
                     0          ,  0, mu,  0,           0,  0, mu,  0,           0,
                     0          ,  0,  0,  0,           0, mu,  0, mu,           0,
                     lambda     ,  0,  0,  0,      lambda,  0,  0,  0, 2*mu+lambda;

  return - dvecF_dx.transpose() * dvecSigma_dvecF * dvecF_dx;
}

Scalar FEM_ElementParameters::compute_volume(const Vec3& x1, const Vec3& x2, const Vec3& x3, const Vec3& x4) const {
  return compute_tetrahedron_volume(x2-x1, x3-x1, x4-x1);
}

void FEM_Element3D::compute_energy_and_derivatives(const PhysicsState& state, EnergyAndDerivatives& out) const {
  // Get the relevant sate
  // ---------------------------------------------------------------
  const Vec3& x1 = p1.get_position(state);
  const Vec3& x2 = p2.get_position(state);
  const Vec3& x3 = p3.get_position(state);
  const Vec3& x4 = p4.get_position(state);

  // Compute deformation, strain and stress tensors
  // ---------------------------------------------------------------
  const Mat3 F = parameters.compute_deformation_tensor(x1, x2, x3, x4);
  const Mat3 epsilon = compute_strain_tensor_from_deformation_tensor(F);
  const Mat3 sigma = parameters.compute_stress_tensor(epsilon);

  // Compute the energy derivatives
  // ---------------------------------------------------------------
  const Scalar volume = abs(parameters.compute_volume(x1, x2, x3, x4));
  const Scalar energy = volume * parameters.get_energy_density(epsilon);
  const Eigen::Vector<Scalar,12> force = volume * parameters.get_force_density(sigma);
  const Eigen::Matrix<Scalar,12,12> df_dx = volume * parameters.get_df_dx_density();

  // Add the energy derivatives to the global structure
  // ---------------------------------------------------------------
  out.energy += energy;

  out.jacobian.segment<3>(p1.index) += force.segment<3>(3*0);
  out.jacobian.segment<3>(p2.index) += force.segment<3>(3*1);
  out.jacobian.segment<3>(p3.index) += force.segment<3>(3*2);
  out.jacobian.segment<3>(p4.index) += force.segment<3>(3*3);

  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      // REVIEW This matrix triplets may be transposed
      // First row
      out.hessian_triplets.emplace_back(p1.index+i, p1.index+j, df_dx(3*0+i, 3*0+j));
      out.hessian_triplets.emplace_back(p1.index+i, p2.index+j, df_dx(3*0+i, 3*1+j));
      out.hessian_triplets.emplace_back(p1.index+i, p3.index+j, df_dx(3*0+i, 3*2+j));
      out.hessian_triplets.emplace_back(p1.index+i, p4.index+j, df_dx(3*0+i, 3*3+j));
      // Second row
      out.hessian_triplets.emplace_back(p2.index+i, p1.index+j, df_dx(3*1+i, 3*0+j));
      out.hessian_triplets.emplace_back(p2.index+i, p2.index+j, df_dx(3*1+i, 3*1+j));
      out.hessian_triplets.emplace_back(p2.index+i, p3.index+j, df_dx(3*1+i, 3*2+j));
      out.hessian_triplets.emplace_back(p2.index+i, p4.index+j, df_dx(3*1+i, 3*3+j));
      // Third row
      out.hessian_triplets.emplace_back(p3.index+i, p1.index+j, df_dx(3*2+i, 3*0+j));
      out.hessian_triplets.emplace_back(p3.index+i, p2.index+j, df_dx(3*2+i, 3*1+j));
      out.hessian_triplets.emplace_back(p3.index+i, p3.index+j, df_dx(3*2+i, 3*2+j));
      out.hessian_triplets.emplace_back(p3.index+i, p4.index+j, df_dx(3*2+i, 3*3+j));
      // Fourth row
      out.hessian_triplets.emplace_back(p4.index+i, p1.index+j, df_dx(3*3+i, 3*0+j));
      out.hessian_triplets.emplace_back(p4.index+i, p2.index+j, df_dx(3*3+i, 3*1+j));
      out.hessian_triplets.emplace_back(p4.index+i, p3.index+j, df_dx(3*3+i, 3*2+j));
      out.hessian_triplets.emplace_back(p4.index+i, p4.index+j, df_dx(3*3+i, 3*3+j));
    }
  }
}

Eigen::Matrix<Scalar,4,3> compute_shape_function_derivative(const Vec3& x1, const Vec3& x2, const Vec3& x3, const Vec3& x4) {
  Eigen::Matrix<Scalar,4,3> ds_dp;
  ds_dp << -1,-1,-1,
            1, 0, 0,
            0, 1, 0,
            0, 0, 1;

  // Here the B matrix is defined as the matrix to transform the coordinate system from the tetrahedron iso-parametric coordinates
  // p' = Bp + x1
  // Where p is a point defined in the iso-paramteric cordinate system and p' is in the same reference frame as x1,x2,etc.
  Mat3 B;
  B << x2-x1, x3-x1, x4-x1;

  return ds_dp * B.inverse();
}

