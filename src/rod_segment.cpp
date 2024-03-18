#include "rod_segment.hpp"
#include <cassert>

inline Eigen::Matrix<Scalar,3,9> dvecR_dtheta_local_matrix(const Mat3& R) {
    return vectorized_levi_civita() * block_matrix(R);
}

Vec3 compute_darboux_vector(const Scalar L0, const Mat3& R1, const Mat3& R2) {
    const Mat3 dR_dx = (R2 - R1) / L0;
    const Mat3 R = (R2 + R1) / 2;
    const Mat3 skew_u = dR_dx * R.transpose();
    const Vec3 u = 0.5 * vectorized_levi_civita() * vectorize_matrix<3>(skew_u);
    return u;
}

Mat3 compute_darboux_vector_theta1_derivative(const Scalar L0, const Mat3& R1, const Mat3& R2) {
    const Mat3 dR_dx = (R2 - R1) / L0;
    const Mat3 R = (R2 + R1) / 2;

    Eigen::Matrix<Scalar,3,9> dvecR_dtheta = 0.5 * dvecR_dtheta_local_matrix(R2);
    const Eigen::Matrix<Scalar,3,9> dvecRx_dtheta = dvecR_dtheta_local_matrix(R2) / L0;

    const Eigen::Matrix<Scalar,9,3> dvec_skewU_dtheta =
        block_matrix<3,3>(R) * dvecRx_dtheta.transpose()                                           // d(dR_dx)/dtheta RT
        + transpose_vectorized_matrix_N<9,3>(block_matrix<3,3>(dR_dx) * dvecR_dtheta.transpose())   // dR/dx (dRT/dtheta)
        ;
    Mat3 du_dtheta = 0.5 * vectorized_levi_civita() * dvec_skewU_dtheta;
    return du_dtheta;
}

Mat3 compute_darboux_vector_theta2_derivative(const Scalar L0, const Mat3& R1, const Mat3& R2) {
    const Mat3 dR_dx = (R2 - R1) / L0;
    const Mat3 R = (R2 + R1) / 2;

    Eigen::Matrix<Scalar,3,9> dvecR_dtheta = 0.5 * dvecR_dtheta_local_matrix(R1);
    const Eigen::Matrix<Scalar,3,9> dvecRx_dtheta = - dvecR_dtheta_local_matrix(R1) / L0;

    const Eigen::Matrix<Scalar,9,3> dvec_skewU_dtheta =
        block_matrix<3,3>(R) * dvecRx_dtheta.transpose()                                           // d(dR_dx)/dtheta RT
        + transpose_vectorized_matrix_N<9,3>(block_matrix<3,3>(dR_dx) * dvecR_dtheta.transpose())   // dR/dx (dRT/dtheta)
        ;
    Mat3 du_dtheta = 0.5 * vectorized_levi_civita() * dvec_skewU_dtheta;
    return du_dtheta;
}

Scalar RodSegmentParameters::get_energy(const Vec3& x1, const Vec3& x2, const Vec3& v1, const Vec3& v2,
                                        const Mat3& R1, const Mat3& R2, const Mat3& R_dot1, const Mat3& R_dot2) const {
    const Scalar one_over_L0 = 1.0 / L0;
    const Vec3 deltaX = x2 - x1;
    const Vec3 u = compute_darboux_vector(L0, R1, R2);

    // Potential energies
    // Stretch
    const Scalar Vs = 0.5 * Ks * std::pow(deltaX.norm() - L0, 2);

    // REVIEW Bending
    const Vec3 deltaU = u - intrinsic_darboux;
    const Scalar Vb = 0.5 * deltaU.transpose() * stiffness_tensor.asDiagonal() * deltaU;

    // Dissipation
    // Translational dissipation
    const Vec3& v_rel = std::pow(one_over_L0, 3) * deltaX * (v2 - v1).dot(deltaX);
    const Scalar Dt = 0.5 * L0 * translational_damping * v_rel.dot(v_rel);

    // TODO Rotational dissipation
    const Scalar Dr = 0.0;

    // Constraint energy
    const Vec3 deltaR_normalized = deltaX.normalized();
    const Mat3 R = 0.5 * (R1 + R2);
    const Vec3 C = (deltaR_normalized - R.col(2)); // We align the 3rd director with the rod segment vector
    const Scalar Ep = 0.5 * L0 * constraint_stiffness * (C).dot(C);

    return Vs + Vb + Dt + Dr;
}

Vec6 RodSegmentParameters::get_energy_gradient(Scalar TimeStep, const Vec3& x1, const Vec3& x2, const Vec3& v1, const Vec3& v2,
                                               const Mat3& R1, const Mat3& R2, const Mat3& R_dot1, const Mat3& R_dot2) const {
    const Scalar one_over_L0 = 1.0 / L0;
    const Scalar one_over_h = 1.0 / TimeStep;
    const Vec3 deltaX = x2 - x1;
    const Scalar L = deltaX.norm();
    const Vec3 u = compute_darboux_vector(L0, R1, R2);

    // Potential energies
    // Stretch
    const Vec3 gradVs = -Ks * (L - L0) * deltaX / L;

    // REVIEW Bending
    const Vec3 deltaU = u - intrinsic_darboux;
    const Vec3 gradVb = deltaU.transpose() * stiffness_tensor.asDiagonal() * compute_darboux_vector_theta1_derivative(L0, R1, R2);

    // Dissipation
    // Translational dissipation
    const Vec3& v_rel = std::pow(one_over_L0, 3) * deltaX * (v2 - v1).dot(deltaX);
    const Vec3 gradDt = L0 * translational_damping * v_rel * one_over_h;

    // TODO Rotational dissipation
    const Vec3 gradDr = Vec3::Zero();

    // Construct the final gradient vector
    Vec6 grad = Vec6::Zero();
    grad << gradVs + gradDt, gradVb + gradDr;
    return grad;
}

Mat6 RodSegmentParameters::get_energy_hessian(Scalar TimeStep, const Vec3& x1, const Vec3& x2, const Vec3& v1, const Vec3& v2,
                                              const Mat3& R1, const Mat3& R2, const Mat3& R_dot1, const Mat3& R_dot2) const {
    const Scalar one_over_L0 = 1.0 / L0;
    const Scalar one_over_h = 1.0 / TimeStep;
    const Vec3 deltaX = x2 - x1;
    const Scalar L = deltaX.norm();
    const Vec3 u = compute_darboux_vector(L0, R1, R2);
    const Mat3 uut = deltaX * deltaX.transpose() / (L*L);

    // Potential energies
    // Stretch
    const Mat3 hessVs = Ks / L * ( (L - L0) * Mat3::Identity() + L0 * uut );

    // REVIEW Bending
    const Vec3 deltaU = u - intrinsic_darboux;
    const Mat3 du_dtheta = compute_darboux_vector_theta1_derivative(L0, R1, R2);
    // Here we are aproximating the hessian!
    const Mat3 hessVb = du_dtheta.transpose() * stiffness_tensor.asDiagonal() * du_dtheta;

    // Dissipation
    // Translational dissipation
    const Mat3 hessDt = L0 * translational_damping * one_over_h * uut;

    // TODO Rotational dissipation
    const Mat3 hessDr = Mat3::Zero();

    // Construct the final matrix:
    Mat6 H = Mat6::Zero();
    H.block<3,3>(0,0) = hessVs + hessDt;
    H.block<3,3>(3,3) = hessVb + hessDr;
    return H;
}


Scalar RodSegment::compute_energy(Scalar TimeStep, const PhysicsState& state) const {
    // Get the relevant sate
    // ---------------------------------------------------------------
    const Vec3 x1 = p1.get_position(state);
    const Vec3 x2 = p2.get_position(state);
    const Vec3 v1 = p1.get_velocity(state);
    const Vec3 v2 = p2.get_velocity(state);
    const Mat3 R1 = rb1.compute_rotation_matrix(state.x);
    const Mat3 R2 = rb2.compute_rotation_matrix(state.x);
    const Mat3 R_dot1 = rb1.compute_rotation_velocity_matrix(TimeStep, state);
    const Mat3 R_dot2 = rb2.compute_rotation_velocity_matrix(TimeStep, state);

    // Compute the energy
    // ---------------------------------------------------------------
    return parameters.get_energy(x1, x2, v1, v2, R1, R2, R_dot1, R_dot2);
}

void RodSegment::compute_energy_and_derivatives(Scalar TimeStep, const PhysicsState& state, EnergyAndDerivatives& out) const {
    // Get the relevant sate
    // ---------------------------------------------------------------
    const Vec3 x1 = p1.get_position(state);
    const Vec3 x2 = p2.get_position(state);
    const Vec3 v1 = p1.get_velocity(state);
    const Vec3 v2 = p2.get_velocity(state);
    const Mat3 R1 = rb1.compute_rotation_matrix(state.x);
    const Mat3 R2 = rb2.compute_rotation_matrix(state.x);
    const Mat3 R_dot1 = rb1.compute_rotation_velocity_matrix(TimeStep, state);
    const Mat3 R_dot2 = rb2.compute_rotation_velocity_matrix(TimeStep, state);

    // Compute the energy derivatives
    // ---------------------------------------------------------------
    const Scalar energy = parameters.get_energy(x1, x2, v1, v2, R1, R2, R_dot1, R_dot2);
    const Vec6 gradient = parameters.get_energy_gradient(TimeStep, x1, x2, v1, v2, R1, R2, R_dot1, R_dot2);
    const Mat6 hessian = parameters.get_energy_hessian(TimeStep, x1, x2, v1, v2, R1, R2, R_dot1, R_dot2);

    // Add the energy derivatives to the global structure
    // ---------------------------------------------------------------

}
