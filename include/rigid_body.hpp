#ifndef RIGID_BODY_H_
#define RIGID_BODY_H_

#include <vector>

#include "linear_algebra.hpp"
#include "physics_state.hpp"

Mat3 compute_initial_inertia_tensor_PARTICLES(Scalar rb_total_mass, const std::vector<Scalar>& vertices);
// Vec3 compute_COM_position(const std::vector<unsigned int>& indices, const std::vector<Scalar>& vertices, RB_MASS_DISTRIBUTION mass_distribution);
Vec3 compute_COM_position_PARTICLES(const std::vector<Scalar>& vertices);
Vec3 compute_COM_position_SHELL(const std::vector<unsigned int>& indices, const std::vector<Scalar>& vertices);
Vec3 compute_COM_position_UNIFORM_VOLUME(const std::vector<unsigned int>& indices, const std::vector<Scalar>& vertices);

struct RigidBody {
    RigidBody(unsigned int index, Scalar mass, Mat3 inertia_tensor0)
        : index(index), mass(mass), inertia_tensor0(inertia_tensor0) {}

    const unsigned int index;
    const Scalar mass;
    const Mat3 inertia_tensor0;

    Vec3 get_COM_position(const PhysicsState& state) const;
    Mat3 compute_rotation_matrix(const Vec& x) const;
    Mat3 compute_inertia_tensor(const Mat3& rotation_matrix) const;
    Mat3 compute_inertia_tensor(const PhysicsState& state) const;
    void update_state(const Vec& dx, Vec& x) const;
};

#endif // RIGID_BODY_H_
