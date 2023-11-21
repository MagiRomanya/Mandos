#ifndef FEM_UNIT_H_
#define FEM_UNIT_H_

#include "linear_algebra.hpp"
#include "particle.hpp"
#include "physics_state.hpp"

struct FEM_ElementParameters {
    FEM_ElementParameters(Scalar mu, Scalar lambda, Vec4 shape_function_coeff, Eigen::Matrix<Scalar,4,3> ds_dx)
    : mu(mu), lambda(lambda), shape_function_coeff(shape_function_coeff), ds_dx(ds_dx) {}

    // Lame coefficients
    const Scalar mu, lambda;

    // Shape function and derivative
    const Vec4 shape_function_coeff;
    const Eigen::Matrix<Scalar,4,3> ds_dx;

    Scalar compute_volume(const Vec3& x1, const Vec3& x2, const Vec3& x3, const Vec3& x4) const;

    Mat3 compute_stress_tensor(const Mat3& sigma) const;

    Scalar get_energy_density(const Mat3& sigma) const;
    Eigen::Vector<Scalar, 12> get_force_density(const Mat3& epsilon) const;
    Eigen::Matrix<Scalar, 12, 12> get_df_dx_density() const;
};

struct FEM_Element3D {
    FEM_Element3D(Particle p1,Particle p2, Particle p3, Particle p4, FEM_ElementParameters param)
        : p1(p1), p2(p2), p3(p3), p4(p4), parameters(param) {}

    const Particle p1, p2, p3, p4;
    const FEM_ElementParameters parameters;
};

#endif // FEM_UNIT_H_
