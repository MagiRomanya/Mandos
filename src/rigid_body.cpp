#include "rigid_body.hpp"
#include "linear_algebra.hpp"

Mat3 skew(const Vec3& v) {
    Mat3 m;
    m <<      0,  -v.z(),   v.y(),
          v.z(),       0,  -v.x(),
         -v.y(),   v.x(),       0;
    return m;
}

Mat3 compute_rotation_matrix_rodrigues(const Vec3& theta) {
    const Scalar angle = theta.norm();
    const Vec3 axis = theta / angle;
    // Rodrigues formula https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula#Matrix_notation
    const Mat3 K = skew(axis);
    return Mat3::Identity() + sin(angle) * K + (1 -cos(angle)) * K * K;
}

Mat3 RigidBody::compute_inertia_tensor(const Mat3& rotation_matrix) const {
    return rotation_matrix * inertia_tensor0 * rotation_matrix.transpose();
}

Scalar RigidBody::get_energy(const Vec3& x, const Vec3& theta) const {
    // TODO
    Scalar energy = 0;
    return energy;
}

Eigen::Vector<Scalar, 6> RigidBody::get_force(const Vec3& x, const Vec3& theta, const Vec3& omega, const Mat3& inertia_tensor) const {
    // Rigid body equations of motion
    // M · a     = F
    // I · alpha = T - omega x I · omega = T_eff

    // TODO handle external forces
    Eigen::Vector<Scalar, 6> f = Eigen::Vector<Scalar, 6>::Zero();
    const Vec3 coriolis_torque = - skew(omega) * inertia_tensor * omega;
    f.tail<3>() += coriolis_torque;

    return f;
}

Eigen::Matrix<Scalar, 6, 6> RigidBody::get_df_dx(const Vec3& x, const Vec3& theta) const {
    // TODO
    Eigen::Matrix<Scalar, 6, 6> df_dx = Eigen::Matrix<Scalar, 6, 6>::Zero();

    return df_dx;
}

void RigidBody::compute_energy_and_derivatives(const PhysicsState& state, EnergyAndDerivatives& out) const {
    // Get the relevant sate
    // ---------------------------------------------------------------
    const Vec3 x = state.x.segment(index, 3);
    const Vec3 v = state.v.segment(index, 3);
    const Vec3 theta = state.x.segment(index + 3, 3);
    const Vec3 omega = state.v.segment(index + 3, 3);

    // Compute rotation matrix and Inerta tensor
    // ---------------------------------------------------------------
    const Mat3 rotation_matrix = compute_rotation_matrix_rodrigues(theta);
    const Mat3 inertia_tensor = compute_inertia_tensor(rotation_matrix);

    // Compute the energy derivatives
    // ---------------------------------------------------------------
    const Scalar energy = get_energy(x, theta);
    const Eigen::Vector<Scalar, 6> force = get_force(x, theta, omega, inertia_tensor);
    const Eigen::Matrix<Scalar, 6, 6> df_dx = get_df_dx(x, theta);

    // Add the energy derivatives to the global structure
    // ---------------------------------------------------------------

    out.energy += energy;
    for (unsigned int i = 0; i<3; i++) {
        out.force[index + i] += force(i);       // force
        out.force[index + 3 + i] += -force(i);  // torque
        // TODO: add force derivatives
        // for (unsigned int j = 0; j<3; j++) {
        //     out.df_dx_triplets.push_back(Triplet(index+i, index+j, df_dx(i, j)));
        // }
    }
}

Vec3 RigidBody::get_COM_position(const PhysicsState& state) const {
    return state.x.segment(index, 3);
}

Mat3 RigidBody::compute_rotation_matrix(const PhysicsState& state) const {
    Vec3 theta = state.x.segment(index+3, 3);
    return compute_rotation_matrix_rodrigues(theta);
}
