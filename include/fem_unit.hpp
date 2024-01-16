#ifndef FEM_UNIT_H_
#define FEM_UNIT_H_

#include "linear_algebra.hpp"
#include "particle.hpp"
#include "physics_state.hpp"

// ds_dx
Eigen::Matrix<Scalar,4,3> compute_shape_function_derivative(const Vec3& x1, const Vec3& x2, const Vec3& x3, const Vec3& x4);

Mat3 compute_deformation_tensor(const Eigen::Matrix<Scalar,9,12>& dvecF_dx, const Vec3& x1, const Vec3& x2, const Vec3& x3, const Vec3& x4);

#define FEM_MATERIAL_MEMBERS \
    MAT(FEM_LinearMaterial, linearMat) \
    MAT(FEM_NeoHookeanMaterial, neoHookMat)

struct FEM_LinearMaterial {
    FEM_LinearMaterial(Scalar mu, Scalar lambda);

    // Lame coefficients
    const Scalar mu, lambda;

    Mat3 compute_strain_tensor(const Mat3& F) const;
    Mat3 compute_stress_tensor(const Mat3& F) const;

    Scalar get_phi(const Mat3& F) const;
    Eigen::Vector<Scalar, 9> get_phi_gradient(const Mat3& sigma) const;
    Eigen::Matrix<Scalar, 9, 9> get_phi_hessian(const Mat3& F) const;
};

struct FEM_NeoHookeanMaterial {
    FEM_NeoHookeanMaterial(Scalar mu, Scalar lambda);

    // Lame coefficients
    const Scalar mu, lambda;

    Mat3 compute_strain_tensor(const Mat3& F) const;
    Mat3 compute_stress_tensor(const Mat3& F) const;

    Scalar get_phi(const Mat3& F) const;
    Eigen::Vector<Scalar, 9> get_phi_gradient(const Mat3& sigma) const;
    Eigen::Matrix<Scalar, 9, 9> get_phi_hessian(const Mat3& F) const;
};

// Eigen::Vector<Scalar, 9> compute_phi_gradinet_finite(const FEM_NeoHookeanMaterial& mat, Scalar dx, const Mat3& F);
// Eigen::Matrix<Scalar, 9, 9> compute_phi_hess_finite(const FEM_NeoHookeanMaterial& mat, Scalar dx, const Mat3& F);

template <typename MaterialType>
struct FEM_Element3D {
    FEM_Element3D(Particle p1,Particle p2, Particle p3, Particle p4, Eigen::Matrix<Scalar, 4, 3> ds_dx, MaterialType material);

    const Particle p1, p2, p3, p4;
    const MaterialType material;

    // Shape derivative (in the form of a block matrix)
    const Eigen::Matrix<Scalar,9,12> dvecF_dx;

    void compute_energy_and_derivatives(Scalar TimeStep, const PhysicsState& state, EnergyAndDerivatives& out) const;

    Scalar compute_volume(const Vec3& x1, const Vec3& x2, const Vec3& x3, const Vec3& x4) const;
};

#endif // FEM_UNIT_H_
