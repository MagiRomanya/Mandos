#include <cmath>
#include <cstring>
#include <ctime>
#include <iostream>
#include <vector>
#include <Eigen/Dense> // for inverse

#include "linear_algebra.hpp"
#include "mandos.hpp"
#include "mesh.hpp"
#include "rigid_body.hpp"
#include "clock.hpp"
#include "utility_functions.hpp"
#include "viewmandos.hpp"

Vec3 rotation_inertia_energy_gradient(const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    const Mat3 J_inertia_tensor = Mat3::Identity();
    const Mat3 R = compute_rotation_matrix_rodrigues(theta);
    const Mat3 R0 = compute_rotation_matrix_rodrigues(theta0);
    const Mat3 Romega = compute_rotation_matrix_rodrigues(TimeStep * omega0);
    // const Mat3 R_guess = (R0 + (R0 - R0old)); // x0 + h* v0
    const Mat3 R_guess = Romega * R0;

    const Mat3 rot_inertia = R * J_inertia_tensor * R_guess.transpose();
    const Mat3 A = (rot_inertia - rot_inertia.transpose()) / 2;
    const Scalar h2 = TimeStep * TimeStep;

    const Vec3 gradient = 2.0f * Vec3(-A(1,2), A(0,2), -A(0,1)) / h2;
    return gradient;
}

Mat3 rotation_inertia_energy_hessian(const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    const Mat3 J_inertia_tensor = Mat3::Identity();
    const Mat3 R = compute_rotation_matrix_rodrigues(theta);
    const Mat3 R0 = compute_rotation_matrix_rodrigues(theta0);
    const Mat3 Romega = compute_rotation_matrix_rodrigues(TimeStep * omega0);
    // const Mat3 R_guess = (R0 + (R0 - R0old)); // x0 + h* v0
    const Mat3 R_guess = Romega * R0;

    const Mat3 rot_inertia = R * J_inertia_tensor * R_guess.transpose();
    const Mat3 S = (rot_inertia + rot_inertia.transpose()) / 2; // Exact hessian
    // const Mat3 S = R * J_inertia_tensor * R.transpose(); // Linear approximation
    const Scalar h2 = TimeStep * TimeStep;

    const Mat3 hessian = 1.0f / h2 * (S.trace() * Mat3::Identity() - S);
    return hessian;
}

Vec3 rotation_inertia_energy_gradient2(const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    const Mat3 J_inertia_tensor = Mat3::Identity();
    const Mat3 R = compute_rotation_matrix_rodrigues(theta);
    const Mat3 R0 = compute_rotation_matrix_rodrigues(theta0);
    const Mat3 Romega = compute_rotation_matrix_rodrigues(TimeStep * omega0);
    // const Mat3 R_guess = (R0 + (R0 - R0old)); // x0 + h* v0
    const Mat3 R_guess = Romega * R0;

    const Mat3 rot_inertia = R * J_inertia_tensor * R_guess.transpose();
    DEBUG_LOG(rot_inertia.transpose());
    std ::cout << "vectorized inertia"
               << " "
               << vectorize_matrix<3>(R_guess).transpose() *
                  block_matrix<3,3>(J_inertia_tensor * R.transpose())
               << std ::endl;
    // const Mat3 A = (rot_inertia - rot_inertia.transpose()) / 2;
    const Mat3 A = (-rot_inertia.transpose()) / 2;
    const Scalar h2 = TimeStep * TimeStep;

    const Vec3 gradient = 1.0f/h2 * vectorized_levi_civita() * vectorize_matrix<3>(A);
    return gradient;
}

Mat3 rotation_inertia_finite_dgradE_dtheta0(const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    Mat3 H;
    const Scalar dx = 0.0001f;
    const Vec3 grad0 = rotation_inertia_energy_gradient2(theta, theta0, omega0, TimeStep);
    for (unsigned i = 0; i < 3; i++) {
        Vec3 dtheta0 = theta0;
        dtheta0[i] += dx;
        const Vec3 grad = rotation_inertia_energy_gradient2(theta, dtheta0, omega0, TimeStep);
        H.row(i) = (grad - grad0) / dx;
    }
    return H;
}

Mat3 rotation_inertia_dgradE_dtheta0(const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    const Mat3 J_inertia_tensor = Mat3::Identity();
    const Mat3 R = compute_rotation_matrix_rodrigues(theta);
    const Mat3 R0 = compute_rotation_matrix_rodrigues(theta0);
    const Mat3 Romega = compute_rotation_matrix_rodrigues(TimeStep * omega0);
    const Mat3 R_guess = Romega * R0;

    const Scalar h2 = TimeStep * TimeStep;
    const Eigen::Matrix<Scalar,3,9> vLeviCivita = vectorized_levi_civita();
    Eigen::Matrix<Scalar,9,3> dvecAdtheta0 = Eigen::Matrix<Scalar,9,3>::Zero();
    // dvecAdtheta0 += -0.5f * block_matrix<3,3>(R*J_inertia_tensor) * (vLeviCivita * block_matrix<3,3>(R_guess.transpose())).transpose();
    dvecAdtheta0 -= 0.5f * block_matrix<3,3>(J_inertia_tensor*R.transpose()) * (vLeviCivita * block_matrix<3,3>(R_guess)).transpose();
    Mat3 H = 1.0f/h2 * vLeviCivita * dvecAdtheta0;
    return H;
}

Eigen::Matrix<Scalar,3,9> dvecR_dtheta_finite(const Vec3& theta) {
    Eigen::Matrix<Scalar,3,9> dvecR_dtheta;
    const Mat3 R0 = compute_rotation_matrix_rodrigues(theta);
    const Eigen::Vector<Scalar,9> vecR0 = vectorize_matrix(R0);
    const Scalar dx = 0.0001f;
    for (unsigned int i = 0; i < 3; i++) {
        Vec3 dtheta = Vec3::Zero();
        dtheta[i] += dx;
        const Eigen::Vector<Scalar,9> vecR = vectorize_matrix<3>(compute_rotation_matrix_rodrigues(dtheta) * R0);
        dvecR_dtheta.row(i) = (vecR - vecR0) / dx;
    }
    return dvecR_dtheta;
}

Eigen::Matrix<Scalar,3,9> dvecR_dtheta_finite2(const Vec3& theta) {
    Eigen::Matrix<Scalar,3,9> dvecR_dtheta;
    const Mat3 R0 = compute_rotation_matrix_rodrigues(theta);
    const Eigen::Vector<Scalar,9> vecR0 = vectorize_matrix(R0);
    const Scalar dx = 0.0001f;
    for (unsigned int i = 0; i < 3; i++) {
        Vec3 dtheta = theta;
        dtheta[i] += dx;
        const Eigen::Vector<Scalar,9> vecR = vectorize_matrix<3>(compute_rotation_matrix_rodrigues(dtheta));
        dvecR_dtheta.row(i) = (vecR - vecR0) / dx;
    }
    return dvecR_dtheta;
}

Eigen::Matrix<Scalar,3,9> dvecR_dtheta_analytic(const Vec3& theta) {
    const Mat3 R = compute_rotation_matrix_rodrigues(theta);
    return vectorized_levi_civita() * block_matrix(R);
}

/*
 * Result = Ra * Rb
*/
inline Vec3 compose_axis_angle(const Vec3& a, const Vec3& b) {
    const Scalar a_angle = a.norm();
    if (a_angle < 1e-8) return b;
    const Vec3 a_axis = a / a_angle;

    const Scalar b_angle = b.norm();
    if (b_angle < 1e-4) return b + a;
    const Vec3 b_axis = b / b_angle;

    // https://math.stackexchange.com/questions/382760/composition-of-two-axis-angle-rotations
    const Scalar sin_a_angle2 = std::sin(a_angle / 2);
    const Scalar cos_a_angle2 = std::cos(a_angle / 2);
    const Scalar sin_b_angle2 = std::sin(b_angle / 2);
    const Scalar cos_b_angle2 = std::cos(b_angle / 2);

    Scalar new_angle = 2.0f * std::acos(cos_a_angle2 * cos_b_angle2
                                      - sin_a_angle2 * sin_b_angle2 * b_axis.dot(a_axis));

    if (fabs(new_angle) < 1e-7) { return Vec3::Zero(); }

    const Vec3 new_axis = 1.0f / std::sin(new_angle/2.0f) * (
        sin_a_angle2 * cos_b_angle2 * a_axis
        + cos_a_angle2 * sin_b_angle2 * b_axis
        + sin_a_angle2 * sin_b_angle2 * cross(a_axis, b_axis));

    new_angle = std::fmod(new_angle, 2.0f * M_PI);
    if (new_angle > M_PI) {
        new_angle -= 2*M_PI;
    }
    return new_angle * new_axis;
}

static const Scalar threshold2 = 1e-3;

inline void compute_axis_angle_jacobian_parts(const Vec3& phi, Mat3& A, Mat3& B) {
    const Mat3 phi_phiT = phi * phi.transpose();
    A = compute_rotation_matrix_rodrigues(phi) - Mat3::Identity() + phi_phiT;
    B = skew(phi) + phi_phiT;
}

inline Mat3 compute_global_to_local_axis_angle_jacobian(const Vec3& phi) {
    if (phi.squaredNorm() < threshold2) return Mat3::Identity();

    Mat3 A, B;
    compute_axis_angle_jacobian_parts(phi, A, B);
    return A.inverse() * B;
}

inline Mat3 compute_local_to_global_axis_angle_jacobian(const Vec3& phi) {
    if (phi.squaredNorm() < threshold2) return Mat3::Identity();

    Mat3 A, B;
    compute_axis_angle_jacobian_parts(phi, A, B);
    return B.inverse() * A;
}

inline Mat3 axis_angle_local_to_global_jacobian_finite(const Vec3& phi0) {
    Mat3 dtheta_dphi;
    const Scalar dx = 0.0001f;
    for (int i = 0; i < 3; i++) {
        Vec3 delta = Vec3::Zero();
        delta[i] = dx;
        const Vec3 dphi = compose_axis_angle(delta, phi0);
        dtheta_dphi.col(i) = (dphi - phi0) / dx;
    }
    return dtheta_dphi;
}

int main(void) {
    const Vec3 phi = M_PI_2 * Vec3(0,1,1).normalized();

    Mat3 finite_jac = axis_angle_local_to_global_jacobian_finite(phi);
    Mat3 a_jac = compute_local_to_global_axis_angle_jacobian(phi);
    Mat3 b_jac = compute_global_to_local_axis_angle_jacobian(phi);
    std ::cout << "finite_jac" << std::endl
               << finite_jac << std ::endl;
    // std ::cout << "a_jac" << std::endl
    //            << a_jac << std ::endl;
    std ::cout << "b_jac" << std::endl
               << b_jac << std ::endl;

    const Scalar angle = phi.norm();
    const Vec3 axis = phi / angle;
    const Scalar half_angle = 0.5 * angle;
    const Mat3 axisaxisT = axis * axis.transpose();
    Mat3 c_jac = half_angle / tan(half_angle) * (Mat3::Identity() - axisaxisT)
        + axisaxisT
        - skew(0.5 * phi);

    std ::cout << "c_jac" << std::endl
               << c_jac << std ::endl;

    Vec3 test = M_PI_2 * Vec3(1,0,1).normalized();
    Mat3 Ra = compute_rotation_matrix_rodrigues(test);
    Vec3 test2 = test + test.normalized() * 2 * M_PI;
    Mat3 Rb = compute_rotation_matrix_rodrigues(test2);
    DEBUG_LOG(Ra);
    DEBUG_LOG(Rb);
    return 0;
}
