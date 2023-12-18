#ifndef RIGID_BODY_H_
#define RIGID_BODY_H_

#include <vector>

#include "linear_algebra.hpp"
#include "physics_state.hpp"

Mat3 compute_initial_inertia_tensor_PARTICLES(Scalar rb_total_mass, const std::vector<Scalar>& vertices);

Vec3 compute_COM_position_PARTICLES(const std::vector<Scalar>& vertices);

Vec3 compute_COM_position_SHELL(const std::vector<unsigned int>& indices, const std::vector<Scalar>& vertices);

Vec3 compute_COM_position_UNIFORM_VOLUME(const std::vector<unsigned int>& indices, const std::vector<Scalar>& vertices);

Vec3 compute_principal_moments_of_inertia(const Mat3& inertia_tensor);

inline Mat3 compute_J_inertia_tensor(const Vec3& I) {
    Mat3 J = Mat3::Zero();
    J(1,1) = - I.x() + I.y() + I.z();
    J(2,2) = + I.x() - I.y() + I.z();
    J(3,3) = + I.x() + I.y() - I.z();
    return 0.5 * J;
}

inline Mat3 compute_J_inertia_tensor(const Mat3& inertia_tensor) {
    const Vec3 I = compute_principal_moments_of_inertia(inertia_tensor);
    return compute_J_inertia_tensor(I);
}

struct RigidBody {
    RigidBody(unsigned int index, Scalar mass, Mat3 inertia_tensor0)
        : index(index), mass(mass), J_inertia_tensor0(inertia_tensor0) {}

    const unsigned int index;
    const Scalar mass;
    const Mat3 J_inertia_tensor0;

    Vec3 get_COM_position(const Vec& x) const;
    Mat3 compute_rotation_matrix(const Vec& x) const;
    Mat3 compute_inertia_tensor(const Mat3& rotation_matrix) const;
    void update_state(const Vec& dx, Vec& x) const;
};

#endif // RIGID_BODY_H_
