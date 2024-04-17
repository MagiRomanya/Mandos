#include <Eigen/Dense> // For inverse matrix

#include "fem_element.hpp"
#include "linear_algebra.hpp"
#include "utility_functions.hpp"

#define MAT(type, name) template struct FEM_Element3D<type>;
FEM_MATERIAL_MEMBERS
#undef MAT

// F == deformation tensor / deformation gradient
// epsilon == strain tensor
// sigma == stress tensor


Mat3 compute_deformation_tensor(const Eigen::Matrix<Scalar,9,12>& dvecF_dx, const Vec3& x1, const Vec3& x2, const Vec3& x3, const Vec3& x4) {
  // X_vec is a 12 component column vector organized as X_vec = x1, x2, x3, x4
  Eigen::Vector<Scalar, 12> X_vec = Eigen::Vector<Scalar, 12>::Zero();
  X_vec.segment<3>(3*0) = x1;
  X_vec.segment<3>(3*1) = x2;
  X_vec.segment<3>(3*2) = x3;
  X_vec.segment<3>(3*3) = x4;
  // F_vec is a column-wise vectorized representation of the F tensor
  // Meaning F_vec = F11, F21, F31, F12, F22, F32, F13, F32, F33
  const Vec9 F_vec = dvecF_dx * X_vec;
  Mat3 F;
  F << F_vec(0), F_vec(3), F_vec(6),
       F_vec(1), F_vec(4), F_vec(7),
       F_vec(2), F_vec(5), F_vec(8);
  return F;
}


FEM_LinearMaterial::FEM_LinearMaterial(Scalar mu, Scalar lambda)
  : mu(mu), lambda(lambda) {}

Mat3 compute_linear_strain_tensor(const Mat3& F) {
  // epsilon
  return 0.5 * (F + F.transpose()) - Mat3::Identity();
}

Scalar FEM_LinearMaterial::get_phi(const Mat3& F) const {
  const Mat3 epsilon = compute_linear_strain_tensor(F);
  return mu * (epsilon.transpose()*epsilon).trace() + 0.5 * lambda * epsilon.trace() * epsilon.trace();
}

Vec9 FEM_LinearMaterial::get_phi_gradient(const Mat3& F) const {
  const Mat3 epsilon = compute_linear_strain_tensor(F);
  // Homogenious material Hook's Law https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
  const Mat3 sigma =  2.0f * mu * epsilon + lambda * epsilon.trace() * Mat3::Identity();
  return vectorize_matrix<3>(sigma);
}

Mat9 FEM_LinearMaterial::get_phi_hessian(const Mat3& F) const {
  Mat9 dvecSigma_dvecF;
  dvecSigma_dvecF << 2*mu+lambda,  0,  0,  0,      lambda,  0,  0,  0,      lambda,
                     0          , mu,  0, mu,           0,  0,  0,  0,           0,
                     0          ,  0, mu,  0,           0,  0, mu,  0,           0,
                     0          , mu,  0, mu,           0,  0,  0,  0,           0,
                     lambda     ,  0,  0,  0, 2*mu+lambda,  0,  0,  0,      lambda,
                     0          ,  0,  0,  0,           0, mu,  0, mu,           0,
                     0          ,  0, mu,  0,           0,  0, mu,  0,           0,
                     0          ,  0,  0,  0,           0, mu,  0, mu,           0,
                     lambda     ,  0,  0,  0,      lambda,  0,  0,  0, 2*mu+lambda;

  return dvecSigma_dvecF;
}

FEM_NeoHookeanMaterial::FEM_NeoHookeanMaterial(Scalar mu_lame, Scalar lambda_lame)
  // mu = 4/3 mu_lame && lambda = lambda_lame + 5/6 mu_lame (eq 16 https://dl.acm.org/doi/pdf/10.1145/3180491)
  : mu(1.33333333333333f * mu_lame), lambda(lambda_lame + 0.833333333333333f * mu_lame)
{
  // Poisson ratio from the neo hookean material parameters
  // const Scalar nu = (lambda - 5.0f / 8.0f *mu) / (2*lambda + mu/4.0f);
}

Scalar FEM_NeoHookeanMaterial::get_phi(const Mat3& F) const {
  const Mat3 C = F.transpose() * F;
  const Scalar IC = C.trace();
  const Scalar J = F.determinant();
  const Scalar alpha = 1.0f + mu / lambda - mu / (4.0f * lambda);
  // const Scalar offset = 9.0f / 32.0f * mu*mu / lambda - mu * std::log(2); // Energy with no deformation
  const Scalar phi = 0.5f * mu * (IC - 3.0f) + 0.5f * lambda * (J - alpha) * (J - alpha) - 0.5f * mu * std::log(IC + 1.0f);
  return phi;
}

template <FEM_Material_Type Material>
Vec9 compute_phi_gradinet_finite(const Material& mat, Scalar dx, const Mat3& F) {
  const Scalar phi0 = mat.get_phi(F);
  Vec9 grad;
  for (unsigned int i = 0; i < 9; i++) {
    Mat3 dF = F;
    const unsigned int row = i % 3;
    const unsigned int col = i / 3;
    dF(row, col) += dx;
    const Scalar dphi = mat.get_phi(dF);
    grad[i] = (dphi - phi0) / dx;
  }
  return grad;
}

template <FEM_Material_Type Material>
Mat9 compute_phi_hess_finite(const Material& mat, Scalar dx, const Mat3& F) {
  const Vec9 grad0 = compute_phi_gradinet_finite(mat, dx, F);
  Mat9 hess;
  for (unsigned int i = 0; i < 9; i++) {
    Mat3 dF = F;
    const unsigned int row = i % 3;
    const unsigned int col = i / 3;
    dF(row, col) += dx;
    const Vec9 dgrad = compute_phi_gradinet_finite(mat, dx, dF);
    hess.col(i) = (dgrad - grad0) / dx;
  }
  return hess;
}

Mat3 determinant_derivative(const Mat3& mat) {
  Mat3 result;
  result << cross(mat.col(1), mat.col(2)), cross(mat.col(2), mat.col(0)), cross(mat.col(0), mat.col(1));
  return result;
}

Vec9 FEM_NeoHookeanMaterial::get_phi_gradient(const Mat3& F) const {
  const Mat3 C = F.transpose() * F;
  const Scalar IC = C.trace();
  const Scalar J = F.determinant();
  const Scalar alpha = 1.0f + mu / lambda - mu / (4.0f * lambda);
  // First Piola-kirchoff Stress (eq 18 https://dl.acm.org/doi/pdf/10.1145/3180491)
  const Mat3 dJ_dF = determinant_derivative(F);
  const Mat3 PK1 = mu * (1.f - 1.f / (IC + 1.f)) * F + lambda * (J - alpha) * dJ_dF;
  return vectorize_matrix<3>(PK1);
}

inline Mat9 compute_vectorized_volume_hessian(const Mat3& F) {
  // (eq 28 https://dl.acm.org/doi/pdf/10.1145/3180491)
  Mat9 vec_H;
  vec_H << Mat3::Zero()    , -skew(F.col(2)) , skew(F.col(1)),
           skew(F.col(2))  , Mat3::Zero()    , -skew(F.col(0)),
           -skew(F.col(1)) , skew(F.col(0))  , Mat3::Zero();
  return vec_H;
}

Mat9 FEM_NeoHookeanMaterial::get_phi_hessian(const Mat3& F) const {
  const Mat3 C = F.transpose() * F;
  const Scalar IC = C.trace();
  const Scalar J = F.determinant();
  const Scalar alpha = 1.0f + mu / lambda - mu / (4.0f * lambda);
  // Energy hessian (eq 22 https://dl.acm.org/doi/pdf/10.1145/3180491)
  // Vectorized energy hessian (eq 41 https://dl.acm.org/doi/pdf/10.1145/3180491)
  const Vec9 vec_F = vectorize_matrix<3>(F);
  const Mat3 dJ_dF = determinant_derivative(F);
  const Vec9 vec_G = vectorize_matrix<3>(dJ_dF);
  const Mat9 vec_H = compute_vectorized_volume_hessian(F);

  Mat9 dvecSigma_dvecF = Mat9::Zero();
  dvecSigma_dvecF += mu * (1.0f - 1.0f / (IC + 1.0f)) * Mat9::Identity();
  dvecSigma_dvecF += mu * 2.0f / ((IC + 1.0f) * (IC + 1.0f)) * vec_F * vec_F.transpose();
  dvecSigma_dvecF += lambda * vec_G * vec_G.transpose();
  dvecSigma_dvecF += lambda * (J - alpha) * vec_H;

// #define ENABLE_PSD_PROJECTION
#ifdef ENABLE_PSD_PROJECTION
  const Eigen::SelfAdjointEigenSolver<Mat9> eigs(dvecSigma_dvecF);
  Vec9 eigenvalues = eigs.eigenvalues();
  Mat9 eigenvectors = eigs.eigenvectors();

  for (int i = 0; i < 9; i++) {
    if (eigenvalues(i) < 0.0) {
      eigenvalues(i) = 0.0;
    }
    }
  dvecSigma_dvecF = eigenvectors * eigenvalues.asDiagonal() * eigenvectors.inverse();
#endif
  return dvecSigma_dvecF;
}

// NOTE We have to transpose the ds_dx matrix before making it a block matrix
template <FEM_Material_Type T>
FEM_Element3D<T>::FEM_Element3D(Particle p1,Particle p2, Particle p3, Particle p4, Eigen::Matrix<Scalar, 4, 3> ds_dx, T material)
  : p1(p1), p2(p2), p3(p3), p4(p4), material(material), dvecF_dx(block_matrix<3,4>(ds_dx.transpose())) {}

template <FEM_Material_Type T>
Scalar FEM_Element3D<T>::compute_volume(const Vec3& x1, const Vec3& x2, const Vec3& x3, const Vec3& x4) const {
  return compute_tetrahedron_volume(x2-x1, x3-x1, x4-x1);
}

template <FEM_Material_Type T>
Scalar FEM_Element3D<T>::compute_energy(Scalar TimeStep, const PhysicsState& state) const {
  // Get the relevant sate
  // ---------------------------------------------------------------
  const Vec3& x1 = p1.get_position(state);
  const Vec3& x2 = p2.get_position(state);
  const Vec3& x3 = p3.get_position(state);
  const Vec3& x4 = p4.get_position(state);

  const Mat3 F = compute_deformation_tensor(dvecF_dx, x1, x2, x3, x4);

  // Compute the energy derivatives
  // ---------------------------------------------------------------
  const Scalar volume = abs(compute_volume(x1, x2, x3, x4));
  const Scalar energy = volume * material.get_phi(F);
  return energy;
}

template <FEM_Material_Type T>
void FEM_Element3D<T>::compute_energy_gradient(Scalar TimeStep, const PhysicsState& state, Vec& grad) const {
  // Get the relevant sate
  // ---------------------------------------------------------------
  const Vec3& x1 = p1.get_position(state);
  const Vec3& x2 = p2.get_position(state);
  const Vec3& x3 = p3.get_position(state);
  const Vec3& x4 = p4.get_position(state);

  // Compute deformation tensor
  // ---------------------------------------------------------------
  const Mat3 F = compute_deformation_tensor(dvecF_dx, x1, x2, x3, x4);

  // Compute the energy derivatives
  // ---------------------------------------------------------------
  const Scalar volume = abs(compute_volume(x1, x2, x3, x4));
  const Eigen::Vector<Scalar,12> gradient = volume * material.get_phi_gradient(F).transpose() *  dvecF_dx;

  grad.segment<3>(p1.index) += gradient.segment<3>(3*0);
  grad.segment<3>(p2.index) += gradient.segment<3>(3*1);
  grad.segment<3>(p3.index) += gradient.segment<3>(3*2);
  grad.segment<3>(p4.index) += gradient.segment<3>(3*3);

}

template <FEM_Material_Type T>
void FEM_Element3D<T>::compute_energy_and_derivatives(Scalar TimeStep, const PhysicsState& state, EnergyAndDerivatives& out) const {
  // Get the relevant sate
  // ---------------------------------------------------------------
  const Vec3& x1 = p1.get_position(state);
  const Vec3& x2 = p2.get_position(state);
  const Vec3& x3 = p3.get_position(state);
  const Vec3& x4 = p4.get_position(state);

  // Compute deformation tensor
  // ---------------------------------------------------------------
  const Mat3 F = compute_deformation_tensor(dvecF_dx, x1, x2, x3, x4);

  // Compute the energy derivatives
  // ---------------------------------------------------------------
  const Scalar volume = abs(compute_volume(x1, x2, x3, x4));
  const Scalar energy = volume * material.get_phi(F);
  const Eigen::Vector<Scalar,12> grad = volume * material.get_phi_gradient(F).transpose() *  dvecF_dx;
  const Eigen::Matrix<Scalar,12,12> hess = volume * dvecF_dx.transpose()* material.get_phi_hessian(F) * dvecF_dx;

  // Add the energy derivatives to the global structure
  // ---------------------------------------------------------------
  out.energy += energy;

  out.gradient.segment<3>(p1.index) += grad.segment<3>(3*0);
  out.gradient.segment<3>(p2.index) += grad.segment<3>(3*1);
  out.gradient.segment<3>(p3.index) += grad.segment<3>(3*2);
  out.gradient.segment<3>(p4.index) += grad.segment<3>(3*3);

  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      // REVIEW This matrix triplets may be transposed
      // -> Not relevant as the hessian should be symmetric
      // First row
      out.hessian_triplets.emplace_back(p1.index+i, p1.index+j, hess(3*0+i, 3*0+j));
      out.hessian_triplets.emplace_back(p1.index+i, p2.index+j, hess(3*0+i, 3*1+j));
      out.hessian_triplets.emplace_back(p1.index+i, p3.index+j, hess(3*0+i, 3*2+j));
      out.hessian_triplets.emplace_back(p1.index+i, p4.index+j, hess(3*0+i, 3*3+j));
      // Second row
      out.hessian_triplets.emplace_back(p2.index+i, p1.index+j, hess(3*1+i, 3*0+j));
      out.hessian_triplets.emplace_back(p2.index+i, p2.index+j, hess(3*1+i, 3*1+j));
      out.hessian_triplets.emplace_back(p2.index+i, p3.index+j, hess(3*1+i, 3*2+j));
      out.hessian_triplets.emplace_back(p2.index+i, p4.index+j, hess(3*1+i, 3*3+j));
      // Third row
      out.hessian_triplets.emplace_back(p3.index+i, p1.index+j, hess(3*2+i, 3*0+j));
      out.hessian_triplets.emplace_back(p3.index+i, p2.index+j, hess(3*2+i, 3*1+j));
      out.hessian_triplets.emplace_back(p3.index+i, p3.index+j, hess(3*2+i, 3*2+j));
      out.hessian_triplets.emplace_back(p3.index+i, p4.index+j, hess(3*2+i, 3*3+j));
      // Fourth row
      out.hessian_triplets.emplace_back(p4.index+i, p1.index+j, hess(3*3+i, 3*0+j));
      out.hessian_triplets.emplace_back(p4.index+i, p2.index+j, hess(3*3+i, 3*1+j));
      out.hessian_triplets.emplace_back(p4.index+i, p3.index+j, hess(3*3+i, 3*2+j));
      out.hessian_triplets.emplace_back(p4.index+i, p4.index+j, hess(3*3+i, 3*3+j));
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

bool is_tetrahedron_inverted(const Vec3& v1, const Vec3& v2, const Vec3& v3, const Vec3& v4) {
  const Vec3 AB = v2 - v1;
  const Vec3 AC = v3 - v1;
  const Vec3 AD = v4 - v1;
  return cross(AB, AC).dot(AD) < 0;
}
