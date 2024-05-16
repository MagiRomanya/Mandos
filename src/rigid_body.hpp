#ifndef RIGID_BODY_H_
#define RIGID_BODY_H_

#include <vector>
#include "linear_algebra.hpp"
#include "physics_state.hpp"
#include "utility_functions.hpp"

/**
 * Clamps the angle part of the axis angle to be between -pi and pi.
 */
Vec3 clamp_axis_angle(const Vec3& axis_angle);

/**
 * Compute the inertia tensor as if the object is made out of mass particles in their vertices.
 *
 * @param rb_total_mass total mass of the object.
 * @param vertices vector of the object's vertices.
 */
Mat3 compute_initial_inertia_tensor_PARTICLES(Scalar rb_total_mass, const std::vector<Scalar>& vertices);

/**
 * Compute the center of mass as if the object is made out of mass particles in their vertices.
 *
 * @param vertices vector of the object's vertices.
 */
Vec3 compute_COM_position_PARTICLES(const std::vector<Scalar>& vertices);

/**
 * Compute the center of mass as if the object was a 2D thin shell (a hollow object)
 *
 * @param indices, vertices vector of the object's indices and vertices.
 */
Vec3 compute_COM_position_SHELL(const std::vector<unsigned int>& indices, const std::vector<Scalar>& vertices);

/**
 * Compute the center of mass as if the object was a 3D volume with uniform volume.
 *
 * @param indices, vertices vector of the object's indices and vertices.
 */
Vec3 compute_COM_position_UNIFORM_VOLUME(const std::vector<unsigned int>& indices, const std::vector<Scalar>& vertices);

/**
 * Computes the pure diagonal components of the inertia tensor.
 */
Vec3 compute_principal_moments_of_inertia(const Mat3& inertia_tensor);

/**
 * Computes the J inertia tensor.
 *
 * @param I principal moments of inertia.
 */
inline Mat3 compute_J_inertia_tensor(const Vec3& I) {
    Mat3 J = Mat3::Zero();
    J(0,0) = - I.x() + I.y() + I.z();
    J(1,1) = + I.x() - I.y() + I.z();
    J(2,2) = + I.x() + I.y() - I.z();
    return 0.5 * J;
}

/**
 * Computes the J inertia tensor.
 *
 * @param inertia_tensor the regular inertia tensor.
 */
inline Mat3 compute_J_inertia_tensor(const Mat3& inertia_tensor) {
    // const Vec3 I = compute_principal_moments_of_inertia(inertia_tensor);
    const Vec3 I = inertia_tensor.diagonal();
    return compute_J_inertia_tensor(I);
}

inline Vec3 compute_principal_moments_of_inertia_from_J(const Mat3& J) {
    return Vec3(J(1,1) + J(2,2), J(2,2) + J(0,0), J(0,0) + J(1,1));
}

Vec3 compose_axis_angle(const Vec3& a, const Vec3& b);

struct RigidBody {
  RigidBody(unsigned int index, Scalar mass, Mat3 inertia_tensor0)
      : index(index), mass(mass),
        J_inertia_tensor0(compute_J_inertia_tensor(inertia_tensor0)) {}

  unsigned int index;
  Scalar mass;
  Mat3 J_inertia_tensor0;

  Vec3 get_COM_position(const Vec &x) const;
  Vec3 get_axis_angle(const Vec &x) const;
  Vec3 compute_angular_momentum(Scalar TimeStep, const PhysicsState &state) const;
  Mat3 compute_rotation_matrix(const Vec &x) const;
  Mat3 compute_rotation_velocity_matrix(const Scalar TimeStep, const PhysicsState &state) const;
};

Mat3 compute_global_to_local_axis_angle_jacobian(const Vec3& phi);

Mat3 compute_local_to_global_axis_angle_jacobian(const Vec3& phi);

/**
 * Compute the derivative of a rotation matrix with respect to local axis angle, evaluated at theta.
 */
inline Eigen::Matrix<Scalar,3,9> dvecR_dtheta_local(const Vec3& theta) {
    const Mat3 R = compute_rotation_matrix_rodrigues(theta);
    return vectorized_levi_civita() * block_matrix<3,3>(R);
}

/**
 * Compute the derivative of a rotation matrix with respect to the global axis angle theta.
 */
inline Eigen::Matrix<Scalar,3,9> dvecR_dtheta_global(const Vec3& theta) {
    const Mat3 jac = compute_local_to_global_axis_angle_jacobian(theta);
    const Mat3 R = compute_rotation_matrix_rodrigues(theta);
    return jac.transpose() * vectorized_levi_civita() * block_matrix<3,3>(R);
}


#endif // RIGID_BODY_H_
