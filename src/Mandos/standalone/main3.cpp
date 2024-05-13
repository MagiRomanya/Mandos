#define _USE_MATH_DEFINES
#include <cmath>
#include <cstring>
#include <ctime>
#include <iostream>
#include <vector>
#include <Eigen/Dense> // for inverse

#include "mandos.hpp"
#include "viewmandos.hpp"
#include "clock.hpp"

const Scalar timestep = 0.1;
const Scalar mass = 1;
const Mat3 inertia_tensor = Mat3::Identity();

void compute_step(const Vec3& x0, const Vec3& v0, Vec3& x1, Vec3& v1) {
    Simulation simulation;
    simulation.TimeStep = timestep;
    RigidBodyHandle rb = RigidBodyHandle(simulation, mass, inertia_tensor)
        .set_initial_orientation(x0)
        .set_initial_angular_velocity(v0)
        ;

    // Compute the simulations tep
    PhysicsState state1 = simulation.initial_state;
    simulation_step(simulation, state1);

    x1 = rb.rb.get_axis_angle(state1.x);
    v1 = rb.rb.get_axis_angle(state1.v);
}

Mat3 compute_dx1_dx0_finite(const Vec3& x0, const Vec3& v0) {
    Vec3 x1, v1;
    compute_step(x0, v0, x1, v1);

    const Scalar dx = 0.00001;
    Mat3 dx1_dx0;
    for (unsigned int i = 0; i < 3; i++) {
        Vec3 dx0 = x0;
        dx0(i) += dx;
        Vec3 dx1, dv1;
        compute_step(dx0, v0, dx1, dv1);
        const Vec3 grad = (dx1 - x1) / dx;
        dx1_dx0.col(i) = grad;
    }
    return dx1_dx0;
}
Mat3 compute_dx1_dv0_finite(const Vec3& x0, const Vec3& v0) {
    Vec3 x1, v1;
    compute_step(x0, v0, x1, v1);

    const Scalar dx = 0.000001;
    Mat3 dx1_dv0;
    for (unsigned int i = 0; i < 3; i++) {
        Vec3 dv0 = v0;
        dv0(i) += dx;
        Vec3 dx1, dv1;
        compute_step(x0, dv0, dx1, dv1);
        const Vec3 grad = (dx1 - x1) / dx;
        dx1_dv0.col(i) = grad;
    }
    return dx1_dv0;
}

Mat3 rotation_inertia_energy_hessian(const Mat3& J_inertia_tensor, const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep);
Vec3 rotation_inertia_energy_gradient2(const Mat3& J_inertia_tensor, const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep);
Mat3 rotation_inertia_finite_dgradE_dtheta(const Mat3& J_inertia_tensor, const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep);
Mat3 rotation_inertia_finite_dgradE_dtheta0(const Mat3& J_inertia_tensor, const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep);
Mat3 rotation_inertia_finite_dgradE_domega0(const Mat3& J_inertia_tensor, const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep);

Mat3 compute_dx1_dx0_ift(const Vec3& x0, const Vec3& v0) {
    Vec3 x1, v1;
    compute_step(x0, v0, x1, v1);

    const Mat3 dgradPsi_dx1 = rotation_inertia_finite_dgradE_dtheta(inertia_tensor, x1, x0, v0, timestep);
    const Mat3 dgradPsi_dx0 = rotation_inertia_finite_dgradE_dtheta0(inertia_tensor, x1, x0, v0, timestep);

    const Mat3 dx1_dx0 = - dgradPsi_dx1.inverse() * dgradPsi_dx0;
    return dx1_dx0;
}
Mat3 compute_dx1_dv0_ift(const Vec3& x0, const Vec3& v0) {
    Vec3 x1, v1;
    compute_step(x0, v0, x1, v1);

    const Mat3 dgradPsi_dx1 = rotation_inertia_finite_dgradE_dtheta(inertia_tensor, x1, x0, v0, timestep);
    const Mat3 dgradPsi_dv0 = rotation_inertia_finite_dgradE_domega0(inertia_tensor, x1, x0, v0, timestep);

    const Mat3 dx1_dv0 = - dgradPsi_dx1.inverse() * dgradPsi_dv0;
    return dx1_dv0;
}

Mat3 compute_dx1_dx0_ift2(const Vec3& x0, const Vec3& v0) {
    Vec3 x1, v1;
    compute_step(x0, v0, x1, v1);

    const Mat3 dgradPsi_dx1 = rotation_inertia_energy_hessian(inertia_tensor, x1, x0, v0, timestep);
    const Mat3 dgradPsi_dx0 = rotation_inertia_finite_dgradE_dtheta0(inertia_tensor, x1, x0, v0, timestep);

    const Mat3 dx1_dx0 = - dgradPsi_dx1.inverse() * dgradPsi_dx0;
    return dx1_dx0;
}

Mat3 compute_dx1_dv0_ift2(const Vec3& x0, const Vec3& v0) {
    Vec3 x1, v1;
    compute_step(x0, v0, x1, v1);

    const Mat3 dgradPsi_dx1 = rotation_inertia_energy_hessian(inertia_tensor, x1, x0, v0, timestep);
    const Mat3 dgradPsi_dv0 = rotation_inertia_finite_dgradE_domega0(inertia_tensor, x1, x0, v0, timestep);

    const Mat3 dx1_dv0 = - dgradPsi_dx1.inverse() * dgradPsi_dv0;
    return dx1_dv0;
}

Scalar loss_function(const Vec3& x) {
    return x.dot(x);
}

Vec3 loss_function_position_derivative(const Vec3& x) {
    return 2 * x;
}

Vec3 compute_loss_position_gradient_finite(const Vec3& x0, const Vec3& v0) {
    Vec3 x1, v1;
    compute_step(x0, v0, x1, v1);
    const Scalar loss0 = loss_function(x1);
    const Scalar dx = 0.000001;
    Vec3 dgdx;
    for (int i = 0; i < 3; i++) {
        Vec3 dx0 = x0;
        dx0(i) += dx;
        Vec3 dx1, dv1;
        compute_step(dx0, v0, dx1, dv1);
        const Scalar loss = loss_function(dx1);
        dgdx(i) = (loss - loss0) / dx;
    }
    return dgdx;
}

int main(void) {
    const Vec3 x0 = Vec3(1,0,0);
    const Vec3 v0 = Vec3(1,2,0);
    Vec3 x1, v1;
    compute_step(x0, v0, x1, v1);
    std ::cout << "x0" << " " << x0.transpose() << std ::endl;
    std ::cout << "x1" << " " << x1.transpose() << std ::endl;

    const Mat3 dx1_dx0_finite = compute_dx1_dx0_finite(x0, v0);
    const Mat3 dx1_dv0_finite = compute_dx1_dv0_finite(x0, v0);

    const Mat3 dx1_dx0_ift = compute_dx1_dx0_ift(x0, v0);
    const Mat3 dx1_dv0_ift = compute_dx1_dv0_ift(x0, v0);

    const Mat3 dx1_dx0_ift2 = compute_dx1_dx0_ift2(x0, v0);
    const Mat3 dx1_dv0_ift2 = compute_dx1_dv0_ift2(x0, v0);

    std ::cout << "dx1_dx0_finite"
               << "\n" << dx1_dx0_finite << std ::endl;
    std ::cout << "dx1_dx0_ift"
               << "\n" << dx1_dx0_ift << std ::endl;
    std ::cout << "dx1_dx0_ift2"
               << "\n" << dx1_dx0_ift2 << std ::endl;
    std ::cout << "dx1_dv0_finite"
               << "\n" << dx1_dv0_finite << std ::endl;
    std ::cout << "dx1_dv0_ift"
               << "\n" << dx1_dv0_ift << std ::endl;
    std ::cout << "dx1_dv0_ift2"
               << "\n" << dx1_dv0_ift2 << std ::endl;

    // Loss gradients
    const Vec3 dgdx_finite = compute_loss_position_gradient_finite(x0, v0);
    const Vec3 partial_dgdx = loss_function_position_derivative(x1);
    const Vec3 dgdx_ift = partial_dgdx.transpose() * dx1_dx0_ift;
    const Vec3 dgdx_ift2 = partial_dgdx.transpose() * dx1_dx0_ift2;
    std ::cout << "partial_dgdx.transpose()"
               << " " << partial_dgdx.transpose() << std ::endl;
    std ::cout << "dgdx_finite"
               << "\n" << dgdx_finite.transpose() << std ::endl;
    std ::cout << "dgdx_ift"
               << "\n" << dgdx_ift.transpose() << std ::endl;
    std ::cout << "dgdx_ift2"
               << "\n" << dgdx_ift2.transpose() << std ::endl;
    return 0;
}


Mat3 rotation_inertia_finite_dgradE_dtheta(const Mat3& J_inertia_tensor, const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    Mat3 H;
    const Scalar dx = 0.0001;
    const Vec3 grad0 = rotation_inertia_energy_gradient2(J_inertia_tensor, theta, theta0, omega0, TimeStep);
    for (unsigned i = 0; i < 3; i++) {
        Vec3 dtheta = theta;
        dtheta[i] += dx;
        const Vec3 grad = rotation_inertia_energy_gradient2(J_inertia_tensor, dtheta, theta0, omega0, TimeStep);
        H.col(i) = (grad - grad0) / dx;
    }
    return H;
}


Mat3 rotation_inertia_finite_dgradE_dtheta0(const Mat3& J_inertia_tensor, const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    Mat3 H;
    const Scalar dx = 0.0001;
    const Vec3 grad0 = rotation_inertia_energy_gradient2(J_inertia_tensor, theta, theta0, omega0, TimeStep);
    for (unsigned i = 0; i < 3; i++) {
        Vec3 dtheta0 = theta0;
        dtheta0[i] += dx;
        const Vec3 grad = rotation_inertia_energy_gradient2(J_inertia_tensor, theta, dtheta0, omega0, TimeStep);
        H.col(i) = (grad - grad0) / dx;
    }
    return H;
}

Mat3 rotation_inertia_finite_dgradE_domega0(const Mat3& J_inertia_tensor, const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    Mat3 H;
    const Scalar dx = 0.0001;
    const Vec3 grad0 = rotation_inertia_energy_gradient2(J_inertia_tensor, theta, theta0, omega0, TimeStep);
    for (unsigned i = 0; i < 3; i++) {
        Vec3 domega0 = omega0;
        domega0[i] += dx;
        const Vec3 grad = rotation_inertia_energy_gradient2(J_inertia_tensor, theta, theta0, domega0, TimeStep);
        H.col(i) = (grad - grad0) / dx;
    }
    return H;
}

Vec3 rotation_inertia_energy_gradient2(const Mat3& J_inertia_tensor, const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    const Mat3 R = compute_rotation_matrix_rodrigues(theta);
    const Mat3 R0 = compute_rotation_matrix_rodrigues(theta0);
    const Mat3 R0old = compute_rotation_matrix_rodrigues(theta0 - omega0 * TimeStep);
    const Mat3 R_guess = (R0 + (R0 - R0old)); // x0 + h* v0

    const Mat3 rot_inertia = R * J_inertia_tensor * R_guess.transpose();
    const Mat3 A = (rot_inertia - rot_inertia.transpose()) / 2;
    const Scalar h2 = TimeStep * TimeStep;

    const Vec3 gradient = 1.0 / h2 * vectorized_levi_civita() * vectorize_matrix<3>(A);
    return gradient;
}

Mat3 rotation_inertia_energy_hessian(const Mat3& J_inertia_tensor, const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    const Mat3 R = compute_rotation_matrix_rodrigues(theta);
    const Mat3 R0 = compute_rotation_matrix_rodrigues(theta0);
    const Mat3 R0old = compute_rotation_matrix_rodrigues(theta0 - omega0 * TimeStep);
    const Mat3 R_guess = (R0 + (R0 - R0old)); // x0 + h* v0

    const Mat3 rot_inertia = R * J_inertia_tensor * R_guess.transpose();
    const Mat3 S = (rot_inertia + rot_inertia.transpose()) / 2; // Exact hessian
    // const Mat3 S = R * J_inertia_tensor * R.transpose(); // Linear approximation
    const Scalar h2 = TimeStep * TimeStep;

    const Mat3 hessian = 1.0 / h2 * (S.trace() * Mat3::Identity() - S);
    return hessian;
}
