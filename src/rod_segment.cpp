#include "rod_segment.hpp"
#include "linear_algebra.hpp"
#include "utility_functions.hpp"
#include <cassert>

inline Eigen::Matrix<Scalar,3,9> dvecR_dtheta_local_matrix(const Mat3& R) {
    return vectorized_levi_civita() * block_matrix<3,3>(R);
}

Vec3 compute_darboux_vector(const Scalar L0, const Mat3& R1, const Mat3& R2) {
    const Mat3 dR_dx = (R2 - R1) / L0;
    const Mat3 R = (R2 + R1) / 2.0;
    const Mat3 skew_u = dR_dx * R.transpose();
    const Vec3 u = 0.5 * vectorized_levi_civita() * vectorize_matrix<3>(skew_u);
    return R.transpose() * u;
    // return u;
}

Mat3 compute_darboux_vector_local_finite_derivative(const Scalar L0, const Mat3& R1, const Mat3& R2) {
    const Vec3 u0 = compute_darboux_vector(L0, R1, R2);
    const Scalar dx = 1e-8;
    Mat3 du_dtheta;
    for (unsigned int i = 0; i < 3; i++) {
        Vec3 theta = Vec3::Zero();
        theta(i) = dx;
        const Mat3 newR1 = compute_rotation_matrix_rodrigues(theta) * R1;
        const Vec3 newU = compute_darboux_vector(L0, newR1, R2);

        du_dtheta.col(i) = (newU - u0) / dx;
    }
    return du_dtheta;
}

Mat3 compute_darboux_vector_local_derivative(const Scalar L0, const Mat3& R1, const Mat3& R2) {
    const Mat3 dR_dx = (R2 - R1) / L0;
    const Mat3 R = (R2 + R1) / 2;

    Eigen::Matrix<Scalar,3,9> dvecR_dtheta = 0.5 * dvecR_dtheta_local_matrix(R1);
    const Eigen::Matrix<Scalar,3,9> dvecRx_dtheta = -dvecR_dtheta_local_matrix(R1) / L0;

    // Compute the derivative of non rotated darboux vector
    const Eigen::Matrix<Scalar,9,3> dvec_skewU_dtheta =
        block_matrix<3,3>(R) * dvecRx_dtheta.transpose()                                            // d(dR_dx)/dtheta RT
        + transpose_vectorized_matrix_N<9,3>(block_matrix<3,3>(dR_dx) * dvecR_dtheta.transpose())   // dR/dx (dRT/dtheta)
        ;
    const Mat3 du_dtheta = 0.5 * vectorized_levi_civita() * dvec_skewU_dtheta;

    // Compute the effect of rotating the darboux with R^T
    const Mat3 skew_u = dR_dx * R.transpose();
    Mat3 dRTu_dtheta = 0.5 * R1.transpose() * skew_u; // Derivative of R^T u wrt theta with u constant

    // Final derivative
    Mat3 result = dRTu_dtheta + R.transpose() * du_dtheta;
    return result;
}

RodSegmentPrecomputedValues::RodSegmentPrecomputedValues(Scalar L0, Scalar TimeStep,
                                                         const Vec3& x1, const Vec3& x2,
                                                         const Vec3& v1, const Vec3& v2,
                                                         const Mat3& R1, const Mat3& R2,
                                                         const Mat3& R_dot1, const Mat3& R_dot2) {
    one_over_L0 = 1.0 / L0;
    one_over_h = 1.0 / TimeStep;
    this->x1 = x1;
    this->x2 = x2;
    this->v1 = v1;
    this->v2 = v2;
    deltaX = x1 - x2;
    L = deltaX.norm();
    one_over_L = 1.0 / L;
    darboux_vector = compute_darboux_vector(L0, R1, R2);
    darboux_vector_derivativeA = compute_darboux_vector_local_derivative(L0, R1, R2);
    darboux_vector_derivativeB = - compute_darboux_vector_local_derivative(L0, R2, R1);
    u = deltaX * one_over_L;
    uut = u * u.transpose();
    v_rel = u * (v1 - v2).dot(u) * one_over_L0;
    this->R1 = R1;
    this->R2 = R2;
    this->R_dot1 = R_dot1;
    this->R_dot2 = R_dot2;
    R = 0.5 * (R1 + R2);
    C = (u - R.col(2));
}

Scalar RodSegmentParameters::compute_energy(const RodSegmentPrecomputedValues& values) const {
    // Potential energies
    // Stretch
    const Scalar Vs = 0.5 * Ks * (values.L - L0) * (values.L - L0);

    // REVIEW Bending
    const Vec3 deltaU = values.darboux_vector - intrinsic_darboux;
    const Scalar Vb = 0.5 * deltaU.transpose() * stiffness_tensor.asDiagonal() * deltaU;

    // Dissipation
    // Translational dissipation
    const Scalar Dt = 0.5 * L0 * translational_damping * values.v_rel.dot(values.v_rel);

    // TODO Rotational dissipation
    const Scalar Dr = 0.0;

    // Constraint energy
    const Scalar Ep = 0.5 * L0 * constraint_stiffness * values.C.squaredNorm();

    return Vs + Vb + Dt + Dr + Ep;
}

Vec3 RodSegmentParameters::compute_energy_linear_gradient(const RodSegmentPrecomputedValues& values) const {
    // Potential energies
    // Stretch
    const Vec3 gradVs = Ks * (values.L - L0) * values.u;

    // Dissipation
    // Translational dissipation
    const Vec3 gradDt = translational_damping * values.uut * (values.v1 - values.v2);

    // Constraint energy
    const Mat3 dL_dx = (Mat3::Identity() - values.uut) * values.one_over_L;
    const Vec3 gradEp_dx = constraint_stiffness * L0 * values.C.transpose() * dL_dx;

    return gradVs + gradDt + gradEp_dx;
}

Vec3 RodSegmentParameters::compute_energy_rotational_gradient_A(const RodSegmentPrecomputedValues& values) const {
    // REVIEW Bending
    const Vec3 deltaU = values.darboux_vector - intrinsic_darboux;
    const Vec3 gradVb = deltaU.transpose() * stiffness_tensor.asDiagonal() * values.darboux_vector_derivativeA;

    // Dissipation
    // TODO Rotational dissipation
    const Vec3 gradDr = Vec3::Zero();

    // Constraint energy
    const Mat3 dd3_dtheta = skew(-0.5 * values.R1.col(2));
    const Vec3 gradEp_dtheta = constraint_stiffness * L0 * values.C.transpose() * (-dd3_dtheta);

    return gradVb + gradDr + gradEp_dtheta;
}

Vec3 RodSegmentParameters::compute_energy_rotational_gradient_B(const RodSegmentPrecomputedValues& values) const {
    // REVIEW Bending
    const Vec3 deltaU = values.darboux_vector - intrinsic_darboux;
    const Mat3 du_dtheta = values.darboux_vector_derivativeB;
    const Vec3 gradVb = deltaU.transpose() * stiffness_tensor.asDiagonal() * du_dtheta;

    // Dissipation
    // TODO Rotational dissipation
    const Vec3 gradDr = Vec3::Zero();

    // Constraint energy
    const Mat3 dd3_dtheta = skew(-0.5 * values.R2.col(2));
    const Vec3 gradEp_dtheta = constraint_stiffness * L0 * values.C.transpose() * (-dd3_dtheta);

    return gradVb + gradDr + gradEp_dtheta;
}

Mat6 RodSegmentParameters::compute_energy_hessian_A(const RodSegmentPrecomputedValues& values) const {
    // Potential energies
    // Stretch
    const Mat3 hessVs = Ks * values.one_over_L * ( (values.L - L0) * Mat3::Identity() + L0 * values.uut );

    // REVIEW Bending
    const Mat3 du_dtheta = values.darboux_vector_derivativeA;
    // Here we are aproximating the hessian!
    const Mat3 hessVb = du_dtheta.transpose() * stiffness_tensor.asDiagonal() * du_dtheta;

    // Dissipation
    // Translational dissipation
    const Mat3 hessDt = translational_damping * values.one_over_h * values.uut;

    // TODO Rotational dissipation
    const Mat3 hessDr = Mat3::Zero();

    // Constraint energy
    const Mat3 dL_dx = (Mat3::Identity() - values.uut) * values.one_over_L;
    const Mat3 dd3_dtheta = - skew(0.5 * values.R1.col(2));
    const Mat3 hessEp_dx2 = constraint_stiffness * L0 * dL_dx.transpose() * dL_dx;
    const Mat3 hessEp_dtheta2 = constraint_stiffness * L0 * dd3_dtheta.transpose() * dd3_dtheta;
    const Mat3 hessEp_dxdtheta = constraint_stiffness * L0 * dL_dx.transpose() * (-dd3_dtheta);


    // Construct the final matrix:
    Mat6 H = Mat6::Zero();
    H.block<3,3>(0, 0) = hessVs  + hessDt + hessEp_dx2;
    H.block<3,3>(0, 3) = hessEp_dxdtheta;
    H.block<3,3>(3, 0) = hessEp_dxdtheta.transpose();
    H.block<3,3>(3, 3) = hessVb + hessDr + hessEp_dtheta2;

    return H;
}

Mat6 RodSegmentParameters::compute_energy_hessian_B(const RodSegmentPrecomputedValues& values) const {
    // Potential energies
    // Stretch
    const Mat3 hessVs = Ks * values.one_over_L * ( (values.L - L0) * Mat3::Identity() + L0 * values.uut );

    // REVIEW Bending
    const Mat3 du_dtheta = values.darboux_vector_derivativeB;
    const Mat3 hessVb = du_dtheta.transpose() * stiffness_tensor.asDiagonal() * du_dtheta;

    // Dissipation
    // Translational dissipation
    const Mat3 hessDt = translational_damping * values.one_over_h * values.uut;

    // TODO Rotational dissipation
    const Mat3 hessDr = Mat3::Zero();

    // Constraint energy
    const Mat3 dL_dx = - (Mat3::Identity() - values.uut) * values.one_over_L;
    const Mat3 dd3_dtheta = skew(-0.5 * values.R2 * Vec3(0,0,1));
    const Mat3 hessEp_dx2 = constraint_stiffness * L0 * dL_dx.transpose() * dL_dx;
    const Mat3 hessEp_dtheta2 = constraint_stiffness * L0 * dd3_dtheta.transpose() * dd3_dtheta;
    const Mat3 hessEp_dxdtheta = constraint_stiffness * L0 * dL_dx.transpose() * (-dd3_dtheta);


    // Construct the final matrix:
    Mat6 H = Mat6::Zero();
    H.block<3,3>(0, 0) = hessVs  + hessDt + hessEp_dx2;
    H.block<3,3>(0, 3) = hessEp_dxdtheta;
    H.block<3,3>(3, 0) = hessEp_dxdtheta.transpose();
    H.block<3,3>(3, 3) = hessVb + hessDr + hessEp_dtheta2;

    return H;
}

Mat6 RodSegmentParameters::compute_energy_hessian_AB(const RodSegmentPrecomputedValues& values) const {
    // Potential energies
    // Stretch
    const Mat3 hessVs = - Ks / values.L * ( (values.L - L0) * Mat3::Identity() + L0 * values.uut );

    // REVIEW Bending
    const Mat3& du_dthetaA = values.darboux_vector_derivativeA;
    const Mat3& du_dthetaB = values.darboux_vector_derivativeB;
    // Here we are aproximating the hessian!
    const Mat3 hessVb = du_dthetaA.transpose() * stiffness_tensor.asDiagonal() * du_dthetaB;

    // Dissipation
    // Translational dissipation
    const Mat3 hessDt = - translational_damping * values.one_over_h * values.uut;

    // TODO Rotational dissipation
    const Mat3 hessDr = Mat3::Zero();

    // Constraint energy
    const Mat3 dL_dxA = (Mat3::Identity() - values.uut) / values.L;
    const Mat3 dL_dxB = - dL_dxA;
    const Mat3 dd3_dthetaA = skew(-0.5 * values.R1.col(2));
    const Mat3 dd3_dthetaB = skew(-0.5 * values.R2.col(2));
    const Mat3 hessEp_dx2 = constraint_stiffness * L0 * dL_dxA.transpose() * dL_dxB;
    const Mat3 hessEp_dtheta2 = constraint_stiffness * L0 * dd3_dthetaA.transpose() * dd3_dthetaB;
    const Mat3 hessEp_dxdtheta = constraint_stiffness * L0 * dL_dxA.transpose() * (-dd3_dthetaB);


    // Construct the final matrix:
    Mat6 H = Mat6::Zero();
    H.block<3,3>(0, 0) = hessVs  + hessDt + hessEp_dx2;
    H.block<3,3>(0, 3) = hessEp_dxdtheta;
    H.block<3,3>(3, 0) = hessEp_dxdtheta;
    H.block<3,3>(3, 3) = hessVb + hessDr + hessEp_dtheta2;
    return H;
}


Scalar RodSegment::compute_energy(Scalar TimeStep, const PhysicsState& state) const {
    // Get the relevant sate
    // ---------------------------------------------------------------
    const Vec3 x1 = rbA.get_COM_position(state.x);
    const Vec3 x2 = rbB.get_COM_position(state.x);
    const Vec3 v1 = rbA.get_COM_position(state.v);
    const Vec3 v2 = rbB.get_COM_position(state.v);
    const Mat3 R1 = rbA.compute_rotation_matrix(state.x);
    const Mat3 R2 = rbB.compute_rotation_matrix(state.x);
    const Mat3 R_dot1 = rbA.compute_rotation_velocity_matrix(TimeStep, state);
    const Mat3 R_dot2 = rbB.compute_rotation_velocity_matrix(TimeStep, state);

    // Compute the energy
    // ---------------------------------------------------------------
    RodSegmentPrecomputedValues values = RodSegmentPrecomputedValues(parameters.L0, TimeStep,
                                                                     x1, x2, v1, v2, R1, R2,
                                                                     R_dot1, R_dot2);
    return parameters.compute_energy(values);
}

void RodSegment::compute_energy_and_derivatives(Scalar TimeStep, const PhysicsState& state, EnergyAndDerivatives& out) const {
    // Get the relevant sate
    // ---------------------------------------------------------------
    const Vec3 x1 = rbA.get_COM_position(state.x);
    const Vec3 x2 = rbB.get_COM_position(state.x);
    const Vec3 v1 = rbA.get_COM_position(state.v);
    const Vec3 v2 = rbB.get_COM_position(state.v);
    const Mat3 R1 = rbA.compute_rotation_matrix(state.x);
    const Mat3 R2 = rbB.compute_rotation_matrix(state.x);
    const Mat3 R_dot1 = rbA.compute_rotation_velocity_matrix(TimeStep, state);
    const Mat3 R_dot2 = rbB.compute_rotation_velocity_matrix(TimeStep, state);

    // Compute the energy derivatives
    // ---------------------------------------------------------------

    RodSegmentPrecomputedValues values = RodSegmentPrecomputedValues(parameters.L0, TimeStep,
                                                                     x1, x2, v1, v2, R1, R2,
                                                                     R_dot1, R_dot2);

    const Scalar energy = parameters.compute_energy(values);
    const Vec3 linear_gradient = parameters.compute_energy_linear_gradient(values);
    const Vec3 rotational_gradient_A = parameters.compute_energy_rotational_gradient_A(values);
    const Vec3 rotational_gradient_B = parameters.compute_energy_rotational_gradient_B(values);

    const Mat6 hessian_A = parameters.compute_energy_hessian_A(values);
    const Mat6 hessian_B = parameters.compute_energy_hessian_B(values);
    const Mat6 hessian_AB = parameters.compute_energy_hessian_AB(values);

    Mat6 hessian_A_finite = Mat6::Zero();
    Mat6 hessian_B_finite = Mat6::Zero();
    Mat6 hessian_AB_finite = Mat6::Zero();
    Mat6 hessian_BA_finite = Mat6::Zero();
    Vec6 grad0A;
    grad0A << linear_gradient, rotational_gradient_A;
    Vec6 grad0B;
    grad0B << -linear_gradient, rotational_gradient_B;
    const Scalar dx = 1e-8;
    for (unsigned int i = 0; i < 3; i++) {
        Vec3 dx_vec = Vec3::Zero();
        dx_vec[i] = dx;
        auto valuesXA = RodSegmentPrecomputedValues(parameters.L0, TimeStep,
                                                    x1 + dx_vec, x2, v1, v2, R1, R2,
                                                    R_dot1, R_dot2);
        auto valuesRA = RodSegmentPrecomputedValues(parameters.L0, TimeStep,
                                                    x1, x2, v1, v2, compute_rotation_matrix_rodrigues(dx_vec)*R1, R2,
                                                    R_dot1, R_dot2);

        auto valuesXB = RodSegmentPrecomputedValues(parameters.L0, TimeStep,
                                                    x1, x2 + dx_vec, v1, v2, R1, R2,
                                                    R_dot1, R_dot2);
        auto valuesRB = RodSegmentPrecomputedValues(parameters.L0, TimeStep,
                                                    x1, x2, v1, v2, R1, compute_rotation_matrix_rodrigues(dx_vec)*R2,
                                                    R_dot1, R_dot2);
        Vec6 dgradXA;
        dgradXA << parameters.compute_energy_linear_gradient(valuesXA), parameters.compute_energy_rotational_gradient_A(valuesXA);
        Vec6 dgradRA;
        dgradRA << parameters.compute_energy_linear_gradient(valuesRA), parameters.compute_energy_rotational_gradient_A(valuesRA);

        Vec6 dgradXB;
        dgradXB << -parameters.compute_energy_linear_gradient(valuesXB), parameters.compute_energy_rotational_gradient_B(valuesXB);
        Vec6 dgradRB;
        dgradRB << -parameters.compute_energy_linear_gradient(valuesRB), parameters.compute_energy_rotational_gradient_B(valuesRB);

        Vec6 dgradXAB;
        dgradXAB << parameters.compute_energy_linear_gradient(valuesXB), parameters.compute_energy_rotational_gradient_A(valuesXB);
        Vec6 dgradRAB;
        dgradRAB << parameters.compute_energy_linear_gradient(valuesRB), parameters.compute_energy_rotational_gradient_A(valuesRB);

        Vec6 dgradXBA;
        dgradXBA << -parameters.compute_energy_linear_gradient(valuesXA), parameters.compute_energy_rotational_gradient_B(valuesXA);
        Vec6 dgradRBA;
        dgradRBA << -parameters.compute_energy_linear_gradient(valuesRA), parameters.compute_energy_rotational_gradient_B(valuesRA);

        hessian_A_finite.col(i) = (dgradXA - grad0A) /dx;
        hessian_A_finite.col(i+3) = (dgradRA - grad0A) / dx;

        hessian_B_finite.col(i) = (dgradXB - grad0B) /dx;
        hessian_B_finite.col(i+3) = (dgradRB - grad0B) / dx;

        hessian_AB_finite.col(i) = (dgradXAB - grad0A) /dx;
        hessian_AB_finite.col(i+3) = (dgradRAB - grad0A) / dx;

        hessian_BA_finite.col(i) = (dgradXBA - grad0B) /dx;
        hessian_BA_finite.col(i+3) = (dgradRBA - grad0B) / dx;
    }
    // std ::cout << "hessian_A"
    //            << "\n" << hessian_A << std ::endl;
    // std ::cout << "hessian_A_finite"
    //            << "\n" << hessian_A_finite << std ::endl;

    // std ::cout << "hessian_B"
    //            << "\n" << hessian_B << std ::endl;
    // std ::cout << "hessian_B_finite"
    //            << "\n" << hessian_B_finite << std ::endl;

    // std ::cout << "hessian_AB"
    //            << "\n" << hessian_AB << std ::endl;
    // std ::cout << "hessian_AB_finite"
    //            << "\n" << hessian_AB_finite << std ::endl;

    // std ::cout << "hessian_AB.transpose()"
    //            << "\n" << hessian_AB.transpose() << std ::endl;
    // std ::cout << "hessian_BA_finite"
    //            << "\n" << hessian_BA_finite << std ::endl;

    // Add the energy derivatives to the global structure
    // ---------------------------------------------------------------
    out.energy += energy;

    // Newton's third law: Every action has an equal and opposite reaction
    out.gradient.segment<3>(rbA.index) += linear_gradient;
    out.gradient.segment<3>(rbB.index) += - linear_gradient;

    // The torque affects differently both rigid bodies
    out.gradient.segment<3>(rbA.index + 3) += rotational_gradient_A;
    out.gradient.segment<3>(rbB.index + 3) += rotational_gradient_B;

    // Fill in the hessian (symetric matrix)
    for (unsigned int i = 0; i < 6; i++) {
        for (unsigned int j = 0; j < 6; j++) {
            out.hessian_triplets.emplace_back(rbA.index + i, rbA.index + j, hessian_A_finite(i,j));
            out.hessian_triplets.emplace_back(rbA.index + i, rbB.index + j, hessian_AB_finite(i,j));
            out.hessian_triplets.emplace_back(rbB.index + i, rbA.index + j, hessian_BA_finite(i,j));
            out.hessian_triplets.emplace_back(rbB.index + i, rbB.index + j, hessian_B_finite(i,j));

            // out.hessian_triplets.emplace_back(rbA.index + i, rbA.index + j, hessian_A(i,j));
            // out.hessian_triplets.emplace_back(rbA.index + i, rbB.index + j, hessian_AB(i,j));
            // out.hessian_triplets.emplace_back(rbB.index + i, rbA.index + j, hessian_AB(j,i));
            // out.hessian_triplets.emplace_back(rbB.index + i, rbB.index + j, hessian_B(i,j));
        }
    }
}
