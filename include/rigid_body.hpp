#ifndef RIGID_BODY_H_
#define RIGID_BODY_H_

#include "linear_algebra.hpp"
#include "physics_state.hpp"

struct RigidBodyParameters {
    Mat3 inertia_tensor0;
};

struct RigidBody {
    RigidBody(unsigned int index, RigidBodyParameters param)
        : index(index), parameters(param) {}

    const unsigned int index;
    const RigidBodyParameters parameters;

    void compute_energy_and_derivatives(const PhysicsState& state, EnergyAndDerivatives& out) const;
    Vec3 get_COM_position(const PhysicsState& state) const;
    Mat3 get_rotation_matrix(const PhysicsState& state) const;

    private:
        Scalar get_energy(const Vec3& x, const Vec3& theta);
        Eigen::Vector<Scalar, 6> get_force(const Vec3& x, const Vec3& theta, const Vec3& omega, const Mat3& inertia_tensor);
        Eigen::Matrix<Scalar, 6, 6> get_df_dx(const Vec3& x, const Vec3& theta);

        Mat3 compute_inertia_tensor(const Mat3& rotation_matrix) const;
};

#endif // RIGID_BODY_H_
