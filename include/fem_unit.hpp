#ifndef FEM_UNIT_H_
#define FEM_UNIT_H_

#include "linear_algebra.hpp"
#include "particle.hpp"
#include "physics_state.hpp"

Eigen::Matrix<Scalar,4,3> compute_shape_function_derivative(const Vec3& x1, const Vec3& x2, const Vec3& x3, const Vec3& x4);

struct FEM_ElementParameters {
    FEM_ElementParameters(Scalar mu, Scalar lambda, Eigen::Matrix<Scalar,4,3> ds_dx);

    // Lame coefficients
    const Scalar mu, lambda;

    // Shape derivative (in the form of a block matrix)
    // const Eigen::Matrix<Scalar,4,3> ds_dx;
    const Eigen::Matrix<Scalar,9,12> dvecF_dx;

    Scalar compute_volume(const Vec3& x1, const Vec3& x2, const Vec3& x3, const Vec3& x4) const;
    Mat3 compute_stress_tensor(const Mat3& sigma) const;
    Mat3 compute_deformation_tensor(const Vec3& x1, const Vec3& x2, const Vec3& x3, const Vec3& x4) const;

    Scalar get_energy_density(const Mat3& sigma) const;
    Eigen::Vector<Scalar, 12> get_force_density(const Mat3& epsilon) const;
    Eigen::Matrix<Scalar, 12, 12> get_df_dx_density() const;
};

struct FEM_Element3D {
    FEM_Element3D(Particle p1,Particle p2, Particle p3, Particle p4, FEM_ElementParameters param)
        : p1(p1), p2(p2), p3(p3), p4(p4), parameters(param) {}

    const Particle p1, p2, p3, p4;
    const FEM_ElementParameters parameters;

    void compute_energy_and_derivatives(const PhysicsState& state, EnergyAndDerivatives& out) const;
};

#endif // FEM_UNIT_H_
