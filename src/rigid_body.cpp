#include <cassert>
#include <vector>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <cmath>

#include "rigid_body.hpp"
#include "linear_algebra.hpp"
#include "sinc.hpp"
#include "utility_functions.hpp"

Mat3 skew(const Vec3& v) {
    Mat3 m;
    m <<     0.,  -v.z(),   v.y(),
          v.z(),      0.,  -v.x(),
         -v.y(),   v.x(),      0.;
    return m;
}

Eigen::Matrix<Scalar,3,4> compute_axis_angle_quaternion_jacobian(const Eigen::Quaternion<Scalar>& q) {
    Eigen::Matrix<Scalar,3,4> dtheta_dq;
    const Scalar q_vec_norm = q.vec().norm();
    const Scalar q_norm2 = q.squaredNorm();
    const Vec3 q_hat = q.vec() / q_vec_norm;
    const Scalar q0 = q.w();
    const Scalar theta = std::atan2(q_vec_norm, q0);
    const Mat3 qqt = q_hat*q_hat.transpose();

    const Mat3 mat = q0 / q_norm2 * qqt + theta * (Mat3::Identity() - qqt) / q_vec_norm;
    const Vec3 vec = - q.vec() / q_norm2;

    dtheta_dq << vec, mat;
    return dtheta_dq;
};

Mat3 compute_rotation_matrix_rodrigues(const Vec3& theta) {
    // https://ipc-sim.github.io/rigid-ipc/assets/rigid_ipc_paper_350ppi.pdf
    // NOTE: There is an error in the paper. The sinc² is multiplied by 2 instead of 1/2.
    const Scalar angle = fmod(theta.norm(), 2 * M_PI);
    const Scalar sinc_angle_2 = sinc(angle / 2.0);
    const Mat3 skew_theta = skew(theta);
    return Mat3::Identity() + sinc(angle) * skew_theta + 0.5 * sinc_angle_2 * sinc_angle_2 * skew_theta * skew_theta;
}

Vec3 RigidBody::get_COM_position(const Vec& x) const {
    return x.segment<3>(index);
}

Vec3 RigidBody::get_axis_angle(const Vec& x) const {
    return x.segment<3>(index + 3);
}

Mat3 RigidBody::compute_rotation_matrix(const Vec& x) const {
    Vec3 theta = get_axis_angle(x);
    return compute_rotation_matrix_rodrigues(theta);
}

Mat3 RigidBody::compute_rotation_velocity_matrix(const Scalar TimeStep, const PhysicsState& state) const {
    const Vec3 theta = get_axis_angle(state.x);
    const Vec3 theta_dot = get_axis_angle(state.v);
    return compute_rotation_matrix_rodrigues(theta - TimeStep * theta_dot) / TimeStep;
}

Scalar compute_tetrahedron_volume(const Vec3& AB, const Vec3& AC, const Vec3& AD) {
    return cross(AB, AC).dot(AD) / 6.;
}

// Computethe volume of the mesh in the same units as the vertex positions.
// This function assumes the mesh to be convex
Scalar compute_mesh_volume(const std::vector<unsigned int>& indices, const std::vector<Scalar>& vertices) {
    Vec3 mesh_center = Vec3::Zero();
    Scalar volume = 0;
    // Compute mesh center
    // ---------------------------------------------------------------
    for (unsigned int v=0; v < vertices.size(); v+=3) {
        mesh_center += Vec3(vertices[v], vertices[v+1], vertices[v+2]);
    }
    mesh_center /= vertices.size() / 3.0;

    // Iterate over triangles and add the volume of each tetrahedron
    // ---------------------------------------------------------------
    for (unsigned int t=0; t < indices.size(); t+=3) {
        const Vec3 v1 = Vec3(vertices[3*indices[t]], vertices[3*indices[t]+1], vertices[3*indices[t]+2]);
        const Vec3 v2 = Vec3(vertices[3*indices[t+1]], vertices[3*indices[t+1]+1], vertices[3*indices[t+1]+2]);
        const Vec3 v3 = Vec3(vertices[3*indices[t+2]], vertices[3*indices[t+2]+1], vertices[3*indices[t+2]+2]);
        volume += abs(compute_tetrahedron_volume(v1 - mesh_center, v2 - mesh_center, v3 - mesh_center));
    }
    return volume;
}

Scalar compute_mesh_surface_area(const std::vector<unsigned int>& indices, const std::vector<Scalar>& vertices) {
    Scalar surface_area = 0.0;
    for (unsigned int t=0; t < indices.size(); t+=3) {
        const Vec3 v1 = Vec3(vertices[3*indices[t]], vertices[3*indices[t]+1], vertices[3*indices[t]+2]);
        const Vec3 v2 = Vec3(vertices[3*indices[t+1]], vertices[3*indices[t+1]+1], vertices[3*indices[t+1]+2]);
        const Vec3 v3 = Vec3(vertices[3*indices[t+2]], vertices[3*indices[t+2]+1], vertices[3*indices[t+2]+2]);
        surface_area += abs(compute_trinagle_area(v2 - v1, v3 - v1));
    }
    return surface_area;
}

Mat3 compute_initial_inertia_tensor_PARTICLES(Scalar rb_total_mass, const std::vector<Scalar>& vertices) {
    Mat3 inertia_tensor = Mat3::Zero();
    const Scalar PARTICLE_MASS = rb_total_mass * 3.0 / vertices.size();
    for (unsigned int i = 0; i < vertices.size(); i+=3) {
        const Vec3 r = Vec3(vertices[i], vertices[i+1], vertices[i+2]);
        const Mat3 skew_r = skew(r);
        inertia_tensor += - PARTICLE_MASS * skew_r * skew_r;
    }
    return inertia_tensor;
}

Vec3 compute_COM_position_PARTICLES(const std::vector<Scalar>& vertices) {
    Vec3 COM_position = Vec3::Zero();
    for (unsigned int i = 0; i < vertices.size(); i+=3) {
        const Vec3 r = Vec3(vertices[i], vertices[i+1], vertices[i+2]);
        COM_position +=r;
    }
    // Average the result
    COM_position *= 3.0 / vertices.size();
    return COM_position;
}

Vec3 compute_COM_position_SHELL(const std::vector<unsigned int>& indices, const std::vector<Scalar>& vertices) {
    Vec3 COM_position = Vec3::Zero();
    Scalar total_surface = 0.0;
    for (unsigned int t=0; t < indices.size(); t+=3) {
        const Vec3 v1 = Vec3(vertices[3*indices[t]], vertices[3*indices[t]+1], vertices[3*indices[t]+2]);
        const Vec3 v2 = Vec3(vertices[3*indices[t+1]], vertices[3*indices[t+1]+1], vertices[3*indices[t+1]+2]);
        const Vec3 v3 = Vec3(vertices[3*indices[t+2]], vertices[3*indices[t+2]+1], vertices[3*indices[t+2]+2]);
        const Scalar triangle_area = abs(compute_trinagle_area(v2 - v1, v3 - v1));
        COM_position += triangle_area * (v1 + v2 + v3) / 3;
        total_surface += total_surface;
    }
    // Average the result
    COM_position /= total_surface;
    return COM_position;
}

Vec3 compute_COM_position_UNIFORM_VOLUME(const std::vector<unsigned int>& indices, const std::vector<Scalar>& vertices) {
    Vec3 COM_position = Vec3::Zero();
    const Vec3 origin = compute_COM_position_PARTICLES(vertices);
    Scalar total_volume = 0.0;
    for (unsigned int t=0; t < indices.size(); t+=3) {
        const Vec3 v1 = Vec3(vertices[3*indices[t]], vertices[3*indices[t]+1], vertices[3*indices[t]+2]);
        const Vec3 v2 = Vec3(vertices[3*indices[t+1]], vertices[3*indices[t+1]+1], vertices[3*indices[t+1]+2]);
        const Vec3 v3 = Vec3(vertices[3*indices[t+2]], vertices[3*indices[t+2]+1], vertices[3*indices[t+2]+2]);
        const Scalar tetrahedron_volume = abs(compute_tetrahedron_volume(v1-origin, v2-origin, v3-origin));
        total_volume += tetrahedron_volume;
        COM_position += tetrahedron_volume * (v1 + v2 + v3 + origin) / 4;
    }
    COM_position /= total_volume;
    return COM_position;
}

/*
 * Result = Ra * Rb
*/
Vec3 compose_axis_angle(const Vec3& a, const Vec3& b) {
    const Scalar a_angle = a.norm();
    if (a_angle < 1e-8) return b;
    const Vec3 a_axis = a / a_angle;

    const Scalar b_angle = b.norm();
    if (b_angle < 1e-8) {
        Scalar a_angle_bounded = std::fmod(a_angle, 2.0 * M_PI);
        if (a_angle_bounded > M_PI) {
            a_angle_bounded -= 2*M_PI;
        }
        return a_angle_bounded * a_axis;
    }
    const Vec3 b_axis = b / b_angle;

    // https://math.stackexchange.com/questions/382760/composition-of-two-axis-angle-rotations
    const Scalar sin_a_angle2 = std::sin(a_angle / 2);
    const Scalar cos_a_angle2 = std::cos(a_angle / 2);
    const Scalar sin_b_angle2 = std::sin(b_angle / 2);
    const Scalar cos_b_angle2 = std::cos(b_angle / 2);

    Scalar new_angle = 2.0 * std::acos(cos_a_angle2 * cos_b_angle2
                                      - sin_a_angle2 * sin_b_angle2 * b_axis.dot(a_axis));

    if (fabs(new_angle) < 1e-8) { return Vec3::Zero(); }

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

inline Mat3 axis_angle_local_to_global_jacobian(const Vec3& theta, const Vec3& dtheta) {
    Mat3 dtheta_dphi;
    const Scalar dx = 0.01;
    const Vec3 phi = compose_axis_angle(dtheta, theta);
    for (int i = 0; i < 3; i++) {
        Vec3 delta = Vec3::Zero();
        delta[i] = dx;
        const Vec3 dphi = compose_axis_angle(dtheta + delta, theta);
        dtheta_dphi.col(i) = (dphi - phi) / dx;
    }
    return dtheta_dphi;
}

Vec3 compute_principal_moments_of_inertia(const Mat3& inertia_tensor) {
    Eigen::EigenSolver<Mat3> solver(inertia_tensor);
    const auto eigenvalues = solver.eigenvalues();
    Mat3 R0 = solver.eigenvectors().real();
    // Ensure that the base rotation is a proper rotation (not a reflection)
    if (R0.determinant() < 0.0) {
        R0.col(0) *= -1.0;
    }
    assert(eigenvalues.imag().isZero());
    return eigenvalues.real();
}

Vec3 clamp_axis_angle(const Vec3& axis_angle) {
    constexpr Scalar tau = 2 * M_PI;
    const Scalar angle = axis_angle.norm();
    if (std::abs(angle) < 0.1) return axis_angle;
    Scalar clamped = std::fmod(angle, tau);
    if (clamped > M_PI) clamped -= tau;
    return axis_angle / angle * clamped;
}

Vec3 RigidBody::compute_angular_momentum(Scalar TimeStep, const PhysicsState& state) const {
    const Vec3 phi = get_axis_angle(state.x);
    const Vec3 phi_dot = get_axis_angle(state.v);
    const Mat3 R = compute_rotation_matrix_rodrigues(phi);
    const Mat3 R0 = compute_rotation_matrix_rodrigues(phi - TimeStep * phi_dot);
    const Mat3 R_dot = (R - R0) / TimeStep;
    const Mat3 W = R.transpose() * R_dot;
    const Mat3 Y = W * J_inertia_tensor0;
    const Mat3 skewY = Y - Y.transpose();
    const Vec3 angular_momentum = Vec3(-skewY(1,2), skewY(0,2), -skewY(0,1));
    return angular_momentum;
}

inline void compute_axis_angle_jacobian_parts(const Vec3& phi, Mat3& A, Mat3& B) {
    const Mat3 phi_phiT = phi * phi.transpose();
    A = compute_rotation_matrix_rodrigues(phi) - Mat3::Identity() + phi_phiT;
    B = skew(phi) + phi_phiT;
}

static const Scalar threshold2 = 1e-5;
Mat3 compute_global_to_local_axis_angle_jacobian(const Vec3& phi) {
    if (phi.squaredNorm() < threshold2) return Mat3::Identity();

    Mat3 A, B;
    compute_axis_angle_jacobian_parts(phi, A, B);
    return A.inverse() * B;
}

Mat3 compute_local_to_global_axis_angle_jacobian(const Vec3& phi) {
    if (phi.squaredNorm() < threshold2) return Mat3::Identity();

    Mat3 A, B;
    compute_axis_angle_jacobian_parts(phi, A, B);
    return B.inverse() * A;
}

