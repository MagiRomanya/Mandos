#include "rigid_body.hpp"
#include "linear_algebra.hpp"

Mat3 skew(const Vec3& v) {
    Mat3 m;
    m <<      0,  -v.z(),   v.y(),
          v.z(),       0,  -v.x(),
         -v.y(),   v.x(),       0;
    return m;
}

Mat3 compute_rotation_matrix(const Vec3& theta) {
    // TODO
    Mat3 rotation_matrix;
    return rotation_matrix;
}

Mat3 RigidBody::compute_inertia_tensor(const Mat3& rotation_matrix) const {
    return rotation_matrix * parameters.inertia_tensor0 * rotation_matrix.transpose();
}

Scalar RigidBody::get_energy(const Vec3& x, const Vec3& theta) {
    // TODO
    Scalar energy;
    return energy;
}

Eigen::Vector<Scalar, 6> RigidBody::get_force(const Vec3& x, const Vec3& theta, const Vec3& omega, const Mat3& inertia_tensor) {
    // TODO
    Eigen::Vector<Scalar, 6> f;

    Vec3 coriolis_torque = - skew(omega) * inertia_tensor * omega;

    return f;
}

Eigen::Matrix<Scalar, 6, 6> RigidBody::get_df_dx(const Vec3& x, const Vec3& theta) {
    // TODO
    Eigen::Vector<Scalar, 6> df_dx;

    return df_dx;
}

void RigidBody::compute_energy_and_derivatives(const PhysicsState& state, EnergyAndDerivatives& out) const {
    // Get the relevant sate
    // ---------------------------------------------------------------
    Vec3 x = state.x.segment(index, 3);
    Vec3 v = state.v.segment(index, 3);
    Vec3 theta = state.x.segment(index + 3, 3);
    Vec3 omega = state.v.segment(index + 3, 3);

    // Compute rotation matrix and Inerta tensor
    // ---------------------------------------------------------------
    Mat3 rotation_matrix = compute_rotation_matrix(theta);
    Mat3 inertia_tensor = compute_inertia_tensor(rotation_matrix);

    // TODO
    // Compute the energy derivatives
    // ---------------------------------------------------------------

    // TODO
    // Add the energy derivatives to the global structure
    // ---------------------------------------------------------------
}

Vec3 RigidBody::get_COM_position(const PhysicsState& state) const {
    return state.x.segment(index, 3);
}

Mat3 RigidBody::get_rotation_matrix(const PhysicsState& state) const {
    Vec3 theta = state.x.segment(index+3, 3);
    return compute_rotation_matrix(theta);
}
