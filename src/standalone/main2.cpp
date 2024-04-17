#include <Eigen/Geometry>
#include "mandos.hpp"
#include "viewmandos.hpp"
#include "../mesh.hpp"
#include "../async_simulation_loop.hpp"


struct RodSegmentStateDEBUG {
    Vec3 x1, x2, v1, v2;
    Mat3 R1, R2, R_dot1, R_dot2;
};

const Vec6 compute_finite_diff_gradientA(const Scalar dx, const RodSegmentParameters& parameters, const RodSegmentStateDEBUG& state) {
    Vec6 grad = Vec6::Zero();
    RodSegmentPrecomputedValues values = RodSegmentPrecomputedValues(parameters.L0, 0.1, state.x1, state.x2, state.v1, state.v2, state.R1, state.R2, state.R_dot1, state.R_dot2);
    const Scalar E0 = parameters.compute_energy(values);
    for (unsigned int i = 0; i < 3; ++i) {
        // Linear grad
        Vec3 dx_vec = Vec3::Zero();
        dx_vec(i) = dx;

        RodSegmentPrecomputedValues valuesX = RodSegmentPrecomputedValues(parameters.L0, 0.1, state.x1 + dx_vec,
                                                                          state.x2, state.v1, state.v2, state.R1, state.R2, state.R_dot1, state.R_dot2);
        const Scalar Ex = parameters.compute_energy(valuesX);
        grad(i) = (Ex - E0) / dx;

        // Rot grad
        const Mat3 newR1 = compute_rotation_matrix_rodrigues(dx_vec) * state.R1;
        RodSegmentPrecomputedValues valuesT = RodSegmentPrecomputedValues(parameters.L0, 0.1, state.x1, state.x2, state.v1, state.v2,
                                                                          newR1,
                                                                          state.R2, state.R_dot1, state.R_dot2);
        const Scalar Etheta = parameters.compute_energy(valuesT);
        grad(i + 3) = (Etheta - E0) / dx;
    }
    return grad;
}

const Vec6 compute_finite_diff_gradientA_axis_angle(const Scalar dx, const RodSegmentParameters& parameters, const RodSegmentStateDEBUG& state) {
    Vec6 grad = Vec6::Zero();
    Eigen::AngleAxis<Scalar> theta1 = Eigen::AngleAxis<Scalar>(state.R1);
    RodSegmentPrecomputedValues values = RodSegmentPrecomputedValues(parameters.L0, 0.1, state.x1, state.x2, state.v1, state.v2,
                                                                     state.R1, state.R2, state.R_dot1, state.R_dot2);
    const Scalar E0 = parameters.compute_energy(values);
    for (unsigned int i = 0; i < 3; ++i) {
        // Linear grad
        Vec3 dx_vec = Vec3::Zero();
        dx_vec(i) = dx;

        RodSegmentPrecomputedValues valuesX = RodSegmentPrecomputedValues(parameters.L0, 0.1, state.x1 + dx_vec,
                                                                          state.x2, state.v1, state.v2, state.R1, state.R2, state.R_dot1, state.R_dot2);
        const Scalar Ex = parameters.compute_energy(valuesX);
        grad(i) = (Ex - E0) / dx;

        // Rot grad
        const Mat3 newR1 = compute_rotation_matrix_rodrigues(theta1.axis()*theta1.angle() + dx_vec);
        RodSegmentPrecomputedValues valuesT = RodSegmentPrecomputedValues(parameters.L0, 0.1, state.x1, state.x2, state.v1, state.v2,
                                                                          newR1,
                                                                          state.R2, state.R_dot1, state.R_dot2);
        const Scalar Etheta = parameters.compute_energy(valuesT);
        grad(i + 3) = (Etheta - E0) / dx;
    }
    return grad;

}

const Vec6 compute_finite_diff_gradientB(const Scalar dx, const RodSegmentParameters& parameters, const RodSegmentStateDEBUG& state) {
    Vec6 grad = Vec6::Zero();
    RodSegmentPrecomputedValues values = RodSegmentPrecomputedValues(parameters.L0, 0.1, state.x1, state.x2, state.v1, state.v2, state.R1, state.R2, state.R_dot1, state.R_dot2);
    const Scalar E0 = parameters.compute_energy(values);
    for (unsigned int i = 0; i < 3; ++i) {
        // Linear grad
        Vec3 dx_vec = Vec3::Zero();
        dx_vec(i) = dx;

        RodSegmentPrecomputedValues valuesX = RodSegmentPrecomputedValues(parameters.L0, 0.1, state.x1,
                                                                          state.x2 + dx_vec,
                                                                          state.v1, state.v2, state.R1, state.R2, state.R_dot1, state.R_dot2);
        const Scalar Ex = parameters.compute_energy(valuesX);
        grad(i) = (Ex - E0) / dx;

        // Rot grad
        const Mat3 newR2 = compute_rotation_matrix_rodrigues(dx_vec) * state.R2;
        RodSegmentPrecomputedValues valuesT = RodSegmentPrecomputedValues(parameters.L0, 0.1, state.x1, state.x2, state.v1, state.v2, state.R1,
                                                                          newR2,
                                                                          state.R_dot1, state.R_dot2);
        const Scalar Etheta = parameters.compute_energy(valuesT);
        grad(i + 3) = (Etheta - E0) / dx;
    }
    return grad;
}

const Mat6 compute_finite_diff_hessianA(const Scalar dx, const RodSegmentParameters& parameters, const RodSegmentStateDEBUG& state) {
    Mat6 hess = Mat6::Zero();
    const Vec6 grad0 = compute_finite_diff_gradientA(dx, parameters, state);
    DEBUG_LOG(grad0.transpose());

    for (unsigned int i = 0; i < 3; ++i) {
        // Linear grad
        Vec3 dx_vec = Vec3::Zero();
        dx_vec(i) = dx;
        RodSegmentStateDEBUG dstate_x = state;
        dstate_x.x1 += dx_vec;
        const Vec6 grad_x = compute_finite_diff_gradientA(dx, parameters, dstate_x);
        hess.col(i) = (grad_x - grad0) / dx;

        // Rot grad
        RodSegmentStateDEBUG dstate_theta = state;
        dstate_theta.R1 = compute_rotation_matrix_rodrigues(dx_vec) * state.R1;
        const Vec6 grad_theta = compute_finite_diff_gradientA(dx, parameters, dstate_theta);
        DEBUG_LOG(grad_theta.transpose());
        hess.col(i + 3) = (grad_theta - grad0) / dx;
    }
    return hess;
}

// Using analytic gradient
const Mat6 compute_finite_diff_hessianA2(const Scalar dx, const RodSegmentParameters& parameters, const RodSegmentStateDEBUG& state) {
    Mat6 hess = Mat6::Zero();
    RodSegmentPrecomputedValues values0(parameters.L0, 0.1, state.x1, state.x2, state.v1, state.v2, state.R1, state.R2, state.R_dot1, state.R_dot2);
    Vec6 grad0;
    Vec3 gradX0 = parameters.compute_energy_linear_gradient(values0);
    Vec3 gradR0 = parameters.compute_energy_rotational_gradient_A(values0);
    grad0 << gradX0, gradR0;

    for (unsigned int i = 0; i < 3; ++i) {
        // Linear grad
        Vec3 dx_vec = Vec3::Zero();
        dx_vec(i) = dx;
        RodSegmentStateDEBUG dstate_x = state;
        dstate_x.x1 += dx_vec;

        RodSegmentPrecomputedValues valuesX(parameters.L0, 0.1, dstate_x.x1, dstate_x.x2, dstate_x.v1, dstate_x.v2, dstate_x.R1, dstate_x.R2, dstate_x.R_dot1, dstate_x.R_dot2);
        Vec6 grad_x;
        Vec3 gradX0 = parameters.compute_energy_linear_gradient(valuesX);
        Vec3 gradR0 = parameters.compute_energy_rotational_gradient_A(valuesX);
        grad_x << gradX0, gradR0;
        hess.col(i) = (grad_x - grad0) / dx;

        // Rot grad
        const Mat3 newR1 = compute_rotation_matrix_rodrigues(dx_vec) * state.R1;
        RodSegmentStateDEBUG dstate_theta = state;
        dstate_theta.R1 = newR1;

        RodSegmentPrecomputedValues valuesR(parameters.L0, 0.1, dstate_theta.x1, dstate_theta.x2, dstate_theta.v1, dstate_theta.v2, dstate_theta.R1, dstate_theta.R2, dstate_theta.R_dot1, dstate_theta.R_dot2);
        Vec6 grad_theta;
        Vec3 gradX1 = parameters.compute_energy_linear_gradient(valuesR);
        Vec3 gradR1 = parameters.compute_energy_rotational_gradient_A(valuesR);
        grad_theta << gradX1, gradR1;
        hess.col(i + 3) = (grad_theta - grad0) / dx;
    }
    return hess;
}

const Mat6 compute_finite_diff_hessianB(const Scalar dx, const RodSegmentParameters& parameters, const RodSegmentStateDEBUG& state) {
    Mat6 hess = Mat6::Zero();
    const Vec6 grad0 = compute_finite_diff_gradientB(dx, parameters, state);

    for (unsigned int i = 0; i < 3; ++i) {
        // Linear grad
        Vec3 dx_vec = Vec3::Zero();
        dx_vec(i) = dx;
        RodSegmentStateDEBUG dstate_x = state;
        dstate_x.x2 += dx_vec;
        const Vec6 grad_x = compute_finite_diff_gradientB(dx, parameters, dstate_x);
        hess.col(i) = (grad_x - grad0) / dx;

        // Rot grad
        const Mat3 newR2 = compute_rotation_matrix_rodrigues(dx_vec) * state.R2;
        RodSegmentStateDEBUG dstate_theta = state;
        dstate_theta.R2 = newR2;
        const Vec6 grad_theta = compute_finite_diff_gradientB(dx, parameters, dstate_theta);
        hess.col(i + 3) = (grad_theta - grad0) / dx;
    }
    return hess;
}

const Mat6 compute_finite_diff_hessianAB(const Scalar dx, const RodSegmentParameters& parameters, const RodSegmentStateDEBUG& state) {
    Mat6 hess = Mat6::Zero();
    const Vec6 grad0 = compute_finite_diff_gradientA(dx, parameters, state);

    for (unsigned int i = 0; i < 3; ++i) {
        // Linear grad
        Vec3 dx_vec = Vec3::Zero();
        dx_vec(i) = dx;
        RodSegmentStateDEBUG dstate_x = state;
        dstate_x.x2 += dx_vec;
        const Vec6 grad_x = compute_finite_diff_gradientA(dx, parameters, dstate_x);
        hess.col(i) = (grad_x - grad0) / dx;

        // Rot grad
        const Mat3 newR2 = compute_rotation_matrix_rodrigues(dx_vec) * state.R2;
        RodSegmentStateDEBUG dstate_theta = state;
        dstate_theta.R2 = newR2;
        const Vec6 grad_theta = compute_finite_diff_gradientA(dx, parameters, dstate_theta);
        hess.col(i + 3) = (grad_theta - grad0) / dx;
    }
    return hess;
}



void dd3_dtheta(const Mat3& R) {
    const Vec3 d3 = R.col(2);
    const Mat3 dd3_dtheta = - skew(d3);
    const Scalar dx = 1e-6;
    Mat3 dd3_dtheta_finite = Mat3::Zero();

    Vec3 grad_finite = Vec3::Zero();
    const Scalar E0 = (Vec3::Ones() - d3).transpose() * (Vec3::Ones() - d3);
    for (int i = 0; i < 3; i++) {
        Vec3 dx_vec = Vec3::Zero();
        dx_vec(i) = dx;
        const Mat3 newR = compute_rotation_matrix_rodrigues(dx_vec) * R;
        dd3_dtheta_finite.col(i) = ( (newR - R) / dx ).col(2);

        const Vec3 dd3 = newR.col(2);
        const Scalar E = (Vec3::Ones() - dd3).transpose() * (Vec3::Ones() - dd3);
        grad_finite(i) = (E - E0) / dx;
    }
    std ::cout << "dd3_dtheta" << "\n" << dd3_dtheta << std ::endl;
    std ::cout << "dd3_dtheta_finite" << "\n" << dd3_dtheta_finite << std ::endl;

    Vec3 grad = - (Vec3::Ones() - d3).transpose() * dd3_dtheta;
    std ::cout << "grad" << "\n" << grad.transpose() << std ::endl;
    std ::cout << "grad_finite" << "\n" << 0.5 * grad_finite.transpose() << std ::endl;
}

void vb_hessian() {
    const Mat3 R1 = compute_rotation_matrix_rodrigues(Vec3(2,6,9).normalized());
    const Mat3 R2 = compute_rotation_matrix_rodrigues(Vec3(2,6,10).normalized());
    // const Mat3 R2 = compute_rotation_matrix_rodrigues(Vec3(-7,-3,4).normalized());

    const Scalar L0 = 1;

    const Vec3 u = compute_darboux_vector(L0, R1, R2);
    const Mat3 du_dtheta = compute_darboux_vector_local_derivative(L0, R1, R2);
    const Vec3 grad0 = u.transpose() * du_dtheta;
    const Scalar dx = 1e-8;

    const Mat3 hess = du_dtheta.transpose() * du_dtheta;
    Mat3 hess_finite = Mat3::Zero();
    Mat3 hess_finite2 = Mat3::Zero();
    Mat3 hess_finite3 = Mat3::Zero();

    for (int i = 0; i < 3; i++) {
        Vec3 dx_vec = Vec3::Zero();
        dx_vec(i) = dx;
        const Mat3 newR1 = compute_rotation_matrix_rodrigues(dx_vec) * R1;
        const Vec3 grad = compute_darboux_vector(L0, newR1, R2).transpose() * du_dtheta;
        hess_finite.col(i) = (grad - grad0) / dx;

        const Vec3 grad2 = u.transpose() * compute_darboux_vector_local_derivative(L0, newR1, R2);
        hess_finite2.col(i) = (grad2 - grad0) / dx;

        const Vec3 grad3 = compute_darboux_vector(L0, newR1, R2).transpose() * compute_darboux_vector_local_derivative(L0, newR1, R2);
        hess_finite3.col(i) = (grad3 - grad0) / dx;

    }
    std ::cout << "hess"
               << "\n" << hess << std ::endl;
    std ::cout << "hess_finite"
               << "\n" << hess_finite << std ::endl;
    std ::cout << "hess_finite2"
               << "\n" << hess_finite2 << std ::endl;
    std ::cout << "hess_finite3"
               << "\n" << hess_finite3 << std ::endl;
}

int main(void) {
    RodSegmentStateDEBUG rod_state = {
    .x1 = Vec3::Zero(),
    .x2 = Vec3(0.0, 0.0, -1.0),
    .v1 = Vec3::Zero(),
    .v2 = Vec3::Zero(),
    .R1 = compute_rotation_matrix_rodrigues(Vec3(1,1,0)),
    .R2 = Mat3::Identity(),
    // .R2 = compute_rotation_matrix_rodrigues(Vec3(0.00,0.2,0)),
    // .R1 = compute_rotation_matrix_rodrigues(Vec3(2,6,9).normalized()),
    // .R2 = compute_rotation_matrix_rodrigues(Vec3(2,6,11).normalized()),
    // .R2 = compute_rotation_matrix_rodrigues(Vec3(-7,-3,4).normalized()),
    .R_dot1 = Mat3::Zero(),
    .R_dot2 = Mat3::Zero(),
    };

    RodSegmentParameters parameters = {
    .Ks = 0.0,
    .L0 = 1.0,
    .translational_damping = 0.0,
    .rotational_damping = 0.0,
    .constraint_stiffness = 100.0,
    .intrinsic_darboux = Vec3::Zero(),
    .stiffness_tensor = 0.0 * Vec3::Ones(),
    };

    // dd3_dtheta(compute_rotation_matrix_rodrigues(Vec3(4, -1, 0.1)));
    // vb_hessian();

    const Scalar dx = 1e-3;
    RodSegmentPrecomputedValues values = RodSegmentPrecomputedValues(parameters.L0, 0.1, rod_state.x1, rod_state.x2, rod_state.v1, rod_state.v2, rod_state.R1, rod_state.R2, rod_state.R_dot1, rod_state.R_dot2);
    const Vec6 grad_finiteA = compute_finite_diff_gradientA(dx, parameters, rod_state);
    const Vec6 grad_finiteA_axis_angle = compute_finite_diff_gradientA_axis_angle(dx, parameters, rod_state);

    const Vec3 gradA_t = parameters.compute_energy_linear_gradient(values);
    const Vec3 gradA_r = parameters.compute_energy_rotational_gradient_A(values);
    Vec6 gradA;
    gradA << gradA_t, gradA_r;
    const Mat6 hess_finiteA = compute_finite_diff_hessianA(dx, parameters, rod_state);
    const Mat6 hessA = parameters.compute_energy_hessian_A(values);
    const Mat6 hess_finiteA2 = compute_finite_diff_hessianA2(dx, parameters, rod_state);

    const Vec6 grad_finiteB = compute_finite_diff_gradientB(dx, parameters, rod_state);
    const Vec3 gradB_t = - parameters.compute_energy_linear_gradient(values);
    const Vec3 gradB_r = parameters.compute_energy_rotational_gradient_B(values);
    Vec6 gradB;
    gradB << gradB_t, gradB_r;
    const Mat6 hess_finiteB = compute_finite_diff_hessianB(dx, parameters, rod_state);
    const Mat6 hessB = parameters.compute_energy_hessian_B(values);

    const Mat6 hess_finiteAB = compute_finite_diff_hessianAB(dx, parameters, rod_state);
    const Mat6 hessAB = parameters.compute_energy_hessian_AB(values);

    std ::cout << "grad_finiteA.transpose()" << "\n" << grad_finiteA.transpose() << std ::endl;
    std ::cout << "grad_finiteA_axis_angle.transpose()" << "\n" << grad_finiteA_axis_angle.transpose() << std ::endl;
    std ::cout << "gradA.transpose()" << "\n" << gradA.transpose() << std ::endl;

    std ::cout << "hess_finiteA" << "\n" << hess_finiteA << std ::endl;
    std ::cout << "hess_finiteA2" << "\n" << hess_finiteA2 << std ::endl;
    std ::cout << "hessA" << "\n" << hessA << std ::endl;

    std ::cout << "grad_finiteB.transpose()" << "\n" << grad_finiteB.transpose() << std ::endl;
    std ::cout << "gradB.transpose()" << "\n" << gradB.transpose() << std ::endl;

    std ::cout << "hess_finiteB" << "\n" << hess_finiteB << std ::endl;
    std ::cout << "hessB" << "\n" << hessB << std ::endl;

    std ::cout << "hess_finiteAB" << "\n" << hess_finiteAB << std ::endl;
    std ::cout << "hessAB" << "\n" << hessAB << std ::endl;
}
