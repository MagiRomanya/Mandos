#include <cassert>
#include <iostream>
#include <vector>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <cmath>

#include "rigid_body.hpp"
#include "linear_algebra.hpp"
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
    Scalar angle = theta.norm();
    angle = std::fmod(angle, 2 * M_PI);

    // Prevent angle = 0, axis not well defined case
    // Using small angle approximation instead
    const Scalar THRESHOLD = 1e-1;
    if (angle < THRESHOLD) return Mat3::Identity() + skew(theta);

    const Vec3 axis = theta / angle;
    // Rodrigues formula https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula#Matrix_notation
    const Mat3 K = skew(axis);
    return Mat3::Identity() + sin(angle) * K + (1. -cos(angle)) * K * K;
}

Mat3 RigidBody::compute_inertia_tensor(const Mat3& rotation_matrix) const {
    return rotation_matrix * J_inertia_tensor0 * rotation_matrix.transpose();
}

Vec3 RigidBody::get_COM_position(const Vec& x) const {
    return x.segment(index, 3);
}

Mat3 RigidBody::compute_rotation_matrix(const Vec& x) const {
    Vec3 theta = x.segment(index+3, 3);
    return compute_rotation_matrix_rodrigues(theta);
}

Scalar compute_tetrahedron_volume(const Vec3& AB, const Vec3& AC, const Vec3& AD) {
    return (skew(AB) * AC).dot(AD) / 6.;
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
        const Mat3 r_star = skew(r);
        inertia_tensor += - PARTICLE_MASS * r_star * r_star;
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

Vec3 update_axis_angle(const Vec3& theta, const Vec3& omega) {
    const Scalar angle = theta.norm();
    if (angle < 1e-4) return theta + omega;
    const Vec3 axis = theta / angle;
    typedef Eigen::Quaternion<Scalar> Quat;
    const Quat q = Quat(Eigen::AngleAxis<Scalar>(angle, axis));
    const Quat q_omega = Quat(0, 0.5 * omega);
    const Quat q_dot = q_omega * q;
    const Quat q_new = Quat(q.w() + q_dot.w(), q.vec() + q_dot.vec());
    const Eigen::AngleAxis<Scalar> angle_axis(q_new);

    const Vec3 new_axis = angle_axis.axis();
    const Scalar new_angle = angle_axis.angle();
    return new_angle * new_axis;
}


void RigidBody::update_state(const Vec& dx, Vec& x) const {
    // Update COM position
    x.segment(index, 3) += dx.segment(index, 3);

    // Update axis-angle rotation
    const Vec3 theta = x.segment(index+3,3);
    const Vec3 omega = dx.segment(index+3,3);
    const Vec3 axis_angle = update_axis_angle(theta, omega);
    x.segment(index+3, 3) = axis_angle;
}

Vec3 compute_principal_moments_of_inertia(const Mat3& inertia_tensor) {
    Eigen::EigenSolver<Mat3> solver(inertia_tensor);
    const auto eigenvalues = solver.eigenvalues();
    assert(eigenvalues.imag().isZero());
    return eigenvalues.real();
}

