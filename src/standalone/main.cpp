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
#include "simulation.hpp"
#include "utility_functions.hpp"
#include "viewmandos.hpp"

Vec3 rotation_inertia_energy_gradient(const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    const Mat3 J_inertia_tensor = Mat3::Identity();
    const Mat3 R = compute_rotation_matrix_rodrigues(theta);
    const Mat3 R0 = compute_rotation_matrix_rodrigues(theta0);
    const Mat3 R0old = compute_rotation_matrix_rodrigues(theta0 - omega0 * TimeStep);
    const Mat3 R_guess = (R0 + (R0 - R0old)); // x0 + h* v0

    const Mat3 rot_inertia = R * J_inertia_tensor * R_guess.transpose();
    const Mat3 A = (rot_inertia - rot_inertia.transpose()) / 2;
    const Scalar h2 = TimeStep * TimeStep;

    const Vec3 gradient = 2.0 * Vec3(-A(1,2), A(0,2), -A(0,1)) / h2;
    return gradient;
}

Mat3 rotation_inertia_energy_hessian(const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    const Mat3 J_inertia_tensor = Mat3::Identity();
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

Vec3 rotation_inertia_energy_gradient2(const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    const Mat3 J_inertia_tensor = Mat3::Identity();
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

Mat3 rotation_inertia_finite_dgradE_dtheta(const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    Mat3 H;
    const Scalar dx = 0.001;
    const Vec3 grad0 = rotation_inertia_energy_gradient(theta, theta0, omega0, TimeStep);
    for (unsigned i = 0; i < 3; i++) {
        Vec3 dtheta = theta0;
        dtheta[i] += dx;
        const Vec3 grad = rotation_inertia_energy_gradient(dtheta, theta0, omega0, TimeStep);
        H.col(i) = (grad - grad0) / dx;
    }
    return H;
}


Mat3 rotation_inertia_finite_dgradE_dtheta0(const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    Mat3 H;
    const Scalar dx = 0.0001;
    const Vec3 grad0 = rotation_inertia_energy_gradient2(theta, theta0, omega0, TimeStep);
    for (unsigned i = 0; i < 3; i++) {
        Vec3 dtheta0 = theta0;
        dtheta0[i] += dx;
        const Vec3 grad = rotation_inertia_energy_gradient2(theta, dtheta0, omega0, TimeStep);
        H.col(i) = (grad - grad0) / dx;
    }
    return H;
}

Mat3 rotation_inertia_finite_dgradE_domega0(const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    Mat3 H;
    const Scalar dx = 0.0001;
    const Vec3 grad0 = rotation_inertia_energy_gradient2(theta, theta0, omega0, TimeStep);
    for (unsigned i = 0; i < 3; i++) {
        Vec3 domega0 = theta0;
        domega0[i] += dx;
        const Vec3 grad = rotation_inertia_energy_gradient2(theta, theta0, domega0, TimeStep);
        H.col(i) = (grad - grad0) / dx;
    }
    return H;
}

Eigen::Matrix<Scalar,3,9> dvecR_dtheta_finite_local(const Vec3& theta) {
    Eigen::Matrix<Scalar,3,9> dvecR_dtheta;
    const Mat3 R0 = compute_rotation_matrix_rodrigues(theta);
    const Eigen::Vector<Scalar,9> vecR0 = vectorize_matrix(R0);
    const Scalar dx = 0.0001;
    for (unsigned int i = 0; i < 3; i++) {
        Vec3 dtheta = Vec3::Zero();
        dtheta[i] += dx;
        const Eigen::Vector<Scalar,9> vecR = vectorize_matrix<3>(compute_rotation_matrix_rodrigues(dtheta) * R0);
        dvecR_dtheta.row(i) = (vecR - vecR0) / dx;
    }
    return dvecR_dtheta;
}

Eigen::Matrix<Scalar,3,9> dvecR_dtheta_finite_global(const Vec3& theta) {
    Eigen::Matrix<Scalar,3,9> dvecR_dtheta;
    const Mat3 R0 = compute_rotation_matrix_rodrigues(theta);
    const Eigen::Vector<Scalar,9> vecR0 = vectorize_matrix(R0);
    const Scalar dx = 0.0001;
    for (unsigned int i = 0; i < 3; i++) {
        Vec3 dtheta = theta;
        dtheta[i] += dx;
        const Eigen::Vector<Scalar,9> vecR = vectorize_matrix<3>(compute_rotation_matrix_rodrigues(dtheta));
        dvecR_dtheta.row(i) = (vecR - vecR0) / dx;
    }
    return dvecR_dtheta;
}

Eigen::Matrix<Scalar,3,9> dvecR_dtheta_analytic_local(const Vec3& theta) {
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

    Scalar new_angle = 2.0 * std::acos(cos_a_angle2 * cos_b_angle2
                                      - sin_a_angle2 * sin_b_angle2 * b_axis.dot(a_axis));

    if (fabs(new_angle) < 1e-7) { return Vec3::Zero(); }

    const Vec3 new_axis = 1.0 / std::sin(new_angle/2.0) * (
        sin_a_angle2 * cos_b_angle2 * a_axis
        + cos_a_angle2 * sin_b_angle2 * b_axis
        + sin_a_angle2 * sin_b_angle2 * cross(a_axis, b_axis));

    new_angle = std::fmod(new_angle, 2.0 * M_PI);
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
    const Scalar dx = 0.0001;
    for (int i = 0; i < 3; i++) {
        Vec3 delta = Vec3::Zero();
        delta[i] = dx;
        const Vec3 dphi = compose_axis_angle(delta, phi0);
        dtheta_dphi.col(i) = (dphi - phi0) / dx;
    }
    return dtheta_dphi;
}

Eigen::Matrix<Scalar,3,9> dvecR_dtheta_analytic_global(const Vec3& theta) {
    // const Mat3 jac = compute_local_to_global_axis_angle_jacobian(theta);
    // return jac.transpose() * dvecR_dtheta_analytic_local(theta);

    const Mat3 jac = compute_local_to_global_axis_angle_jacobian(theta);
    const Mat3 R = compute_rotation_matrix_rodrigues(theta);
    return jac.transpose() * vectorized_levi_civita() * block_matrix(R);
}

Mat3 rotation_inertia_dgradE_dtheta0(const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    const Mat3 J_inertia_tensor = Mat3::Identity();
    const Mat3 R = compute_rotation_matrix_rodrigues(theta);

    const Scalar h2 = TimeStep * TimeStep;
    const Eigen::Matrix<Scalar,3,9> vLeviCivita = vectorized_levi_civita();
    const Eigen::Matrix<Scalar,3,9> dvecRguess_dtheta0 = 2 * dvecR_dtheta_analytic_global(theta0) - dvecR_dtheta_analytic_global(theta0 - omega0 * TimeStep);

    const Eigen::Matrix<Scalar,9,3> dvecRMR_guess_dtheta0 = block_matrix<3,3>(R * J_inertia_tensor) * dvecRguess_dtheta0.transpose();;
    const Eigen::Matrix<Scalar,9,3> dvecAdtheta0 = 0.5 * (transpose_vectorized_matrix(dvecRMR_guess_dtheta0) - dvecRMR_guess_dtheta0);

    Mat3 H = 1.0 / h2 * vLeviCivita * dvecAdtheta0;
    return H;
}
Mat3 rotation_inertia_dgradE_domega0(const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    const Mat3 J_inertia_tensor = Mat3::Identity();
    const Mat3 R = compute_rotation_matrix_rodrigues(theta);

    const Scalar h2 = TimeStep * TimeStep;
    const Eigen::Matrix<Scalar,3,9> vLeviCivita = vectorized_levi_civita();
    const Eigen::Matrix<Scalar,3,9> dvecRguess_domega0 = TimeStep * dvecR_dtheta_analytic_global(theta0 - omega0 * TimeStep);

    const Eigen::Matrix<Scalar,9,3> dvecRMR_guess_domega0 = block_matrix<3,3>(R * J_inertia_tensor) * dvecRguess_domega0.transpose();;
    const Eigen::Matrix<Scalar,9,3> dvecAdtheta0 = 0.5 * (transpose_vectorized_matrix(dvecRMR_guess_domega0) - dvecRMR_guess_domega0);

    Mat3 H = 1.0 / h2 * vLeviCivita * dvecAdtheta0;
    return H;
}
Mat3 rotation_inertia_dgradE_dtheta(const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    const Mat3 J_inertia_tensor = Mat3::Identity();
    const Mat3 R0 = compute_rotation_matrix_rodrigues(theta0);
    const Mat3 R0old = compute_rotation_matrix_rodrigues(theta0 - TimeStep * omega0);
    const Mat3 Rguess = (R0 + (R0 - R0old)); // x0 + h* v0

    const Scalar h2 = TimeStep * TimeStep;
    const Eigen::Matrix<Scalar,3,9> vLeviCivita = vectorized_levi_civita();
    const Eigen::Matrix<Scalar,3,9> dvecR_dtheta = dvecR_dtheta_analytic_global(theta);

    const Eigen::Matrix<Scalar,9,3> dvecRMR_guess_dtheta = block_matrix<3,3>(Rguess * J_inertia_tensor) * dvecR_dtheta.transpose();;
    const Eigen::Matrix<Scalar,9,3> dvecAdtheta = 0.5 * (dvecRMR_guess_dtheta - transpose_vectorized_matrix(dvecRMR_guess_dtheta));

    Mat3 H = 1.0 / h2 * vLeviCivita * dvecAdtheta;
    return H;
}



inline Mat3 compute_R_guess(const Mat3& R0, const Mat3 R0old) {
    return (R0 + (R0 - R0old)); // x0 + h* v0
}

Eigen::Matrix<Scalar,3,9> compute_dR_guess_dtheta0_finite(const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    const Scalar dx = 0.000000001;
    const Mat3 R0 = compute_rotation_matrix_rodrigues(theta0);
    const Mat3 R0old = compute_rotation_matrix_rodrigues(theta0 - omega0 * TimeStep);
    const Mat3 R_guess = compute_R_guess(R0, R0old);
    const Vec9 vR_guess0 = vectorize_matrix(R_guess);
    Eigen::Matrix<Scalar,3,9> dR_dtheta0;
    for (int i = 0; i < 3; i++) {
        Vec3 dtheta0 = theta0;
        dtheta0(i) += dx;
        const Mat3 R0 = compute_rotation_matrix_rodrigues(dtheta0);
        const Mat3 R0old = compute_rotation_matrix_rodrigues(dtheta0 - omega0 * TimeStep);
        const Mat3 R_guess = compute_R_guess(R0, R0old);
        const Vec9 vR_guess = vectorize_matrix(R_guess);
        dR_dtheta0.row(i) = (vR_guess - vR_guess0) / dx;
    }
    return dR_dtheta0;
}

Eigen::Matrix<Scalar,3,9> compute_dR_guessT_dtheta0_finite(const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    const Scalar dx = 0.000000001;
    const Mat3 R0 = compute_rotation_matrix_rodrigues(theta0);
    const Mat3 R0old = compute_rotation_matrix_rodrigues(theta0 - omega0 * TimeStep);
    const Mat3 R_guess = compute_R_guess(R0, R0old).transpose();
    const Vec9 vR_guess0 = vectorize_matrix(R_guess);
    Eigen::Matrix<Scalar,3,9> dR_dtheta0;
    for (int i = 0; i < 3; i++) {
        Vec3 dtheta0 = theta0;
        dtheta0(i) += dx;
        const Mat3 R0 = compute_rotation_matrix_rodrigues(dtheta0);
        const Mat3 R0old = compute_rotation_matrix_rodrigues(dtheta0 - omega0 * TimeStep);
        const Mat3 R_guess = compute_R_guess(R0, R0old).transpose();
        const Vec9 vR_guess = vectorize_matrix(R_guess);
        dR_dtheta0.row(i) = (vR_guess - vR_guess0) / dx;
    }
    return dR_dtheta0;
}

Eigen::Matrix<Scalar,3,9> compute_dR_guess_dtheta0(const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    return 2 * dvecR_dtheta_analytic_global(theta0) - dvecR_dtheta_analytic_global(theta0 - omega0 * TimeStep);
}


int main(void) {
    const Vec3 theta = M_PI_2 * Vec3(0,0,1.1).normalized();
    const Vec3 theta0 = M_PI_2 * Vec3(0,0,1).normalized();
    const Vec3 omega0 = M_PI_2 * Vec3(0,0,1).normalized();
    const Scalar TimeStep = 0.1;

    const Mat3 dgradE_dtheta = rotation_inertia_dgradE_dtheta(theta, theta0, omega0, TimeStep);
    const Mat3 dgradE_dtheta_finite = rotation_inertia_finite_dgradE_dtheta(theta, theta0, omega0, TimeStep);
    const Mat3 dgradE_dtheta0 = rotation_inertia_dgradE_dtheta0(theta, theta0, omega0, TimeStep);
    const Mat3 dgradE_domega0 = rotation_inertia_dgradE_domega0(theta, theta0, omega0, TimeStep);
    const Mat3 dgradE_dtheta0_finite = rotation_inertia_finite_dgradE_dtheta0(theta, theta0, omega0, TimeStep);
    const Mat3 dgradE_domega0_finite = rotation_inertia_finite_dgradE_domega0(theta, theta0, omega0, TimeStep);
    const Mat3 jac0 = compute_global_to_local_axis_angle_jacobian(theta0).inverse();
    const Mat3 jac = compute_local_to_global_axis_angle_jacobian(theta);

    std ::cout << "Hess"
               << "\n" << rotation_inertia_energy_hessian(theta, theta0, omega0, TimeStep) << std ::endl;
    std ::cout << "Rot Hess"
               << "\n" << jac.transpose() * rotation_inertia_energy_hessian(theta, theta0, omega0, TimeStep) << std ::endl;
    std ::cout << "dgradE_dtheta_finite"
               << "\n" << dgradE_dtheta_finite << std ::endl;
    std ::cout << "dgradE_dtheta"
               << "\n" << dgradE_dtheta << std ::endl;
    std ::cout << "dgradE_dtheta0_finite"
               << "\n" << dgradE_dtheta0_finite << std ::endl;
    std ::cout << "dgradE_dtheta0"
               << "\n" << dgradE_dtheta0 << std ::endl;
    std ::cout << "dgradE_domega0"
               << "\n" << dgradE_domega0 << std ::endl;
    std ::cout << "dgradE_domega0_finite"
               << "\n" << dgradE_domega0_finite << std ::endl;

    DEBUG_LOG(rotation_inertia_energy_gradient(theta, theta0, theta0, TimeStep).transpose());
    DEBUG_LOG(rotation_inertia_energy_gradient2(theta, theta0, theta0, TimeStep).transpose());

    std ::cout << "dvecR_dtheta_finite_local(theta)"
               << "\n" << dvecR_dtheta_finite_local(theta) << std ::endl;
    std ::cout << "dvecR_dtheta_analytic(theta)"
               << "\n" << dvecR_dtheta_analytic_local(theta) << std ::endl;
    std ::cout << "dvecR_dtheta_analytic_global(theta)"
               << "\n" << dvecR_dtheta_analytic_global(theta) << std ::endl;
    std ::cout << "dvecR_dtheta_finite_global(theta)"
               << "\n" << dvecR_dtheta_finite_global(theta) << std ::endl;

    std ::cout << "compute_dR_guess_dtheta0_finite(theta0, omega0, TimeStep)"
               << "\n"
               << compute_dR_guess_dtheta0_finite(theta0, omega0, TimeStep)
               << std ::endl;
    std ::cout << "compute_dR_guess_dtheta0(theta0, omega0, TimeStep)"
               << "\n" << compute_dR_guess_dtheta0(theta0, omega0, TimeStep)
               << std ::endl;

    std ::cout << "compute_dRT_guess_dtheta0(theta0, omega0, TimeStep)"
    << "\n" << transpose_vectorized_matrix(compute_dR_guess_dtheta0(theta0, omega0, TimeStep))
               << std ::endl;
    std ::cout << "compute_dR_guessT_dtheta0_finite(theta0, omega0, TimeStep)"
               << "\n"
               << compute_dR_guessT_dtheta0_finite(theta0, omega0, TimeStep)
               << std ::endl;

    return 0;
}
