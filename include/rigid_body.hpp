#ifndef RIGID_BODY_H_
#define RIGID_BODY_H_

#include "linear_algebra.hpp"
#include "physics_state.hpp"

struct RigidBody {
    RigidBody(unsigned int index, Scalar mass, Mat3 inertia_tensor0)
        : index(index), inertia_tensor0(inertia_tensor0), mass(mass) {}

    const unsigned int index;
    const Scalar mass;
    const Mat3 inertia_tensor0;

    void compute_energy_and_derivatives(const PhysicsState& state, EnergyAndDerivatives& out) const;
    Vec3 get_COM_position(const PhysicsState& state) const;
    Mat3 compute_rotation_matrix(const PhysicsState& state) const;
    Mat3 compute_inertia_tensor(const Mat3& rotation_matrix) const;

    private:
        Scalar get_energy(const Vec3& x, const Vec3& theta) const;
        Eigen::Vector<Scalar, 6> get_force(const Vec3& x, const Vec3& theta, const Vec3& omega, const Mat3& inertia_tensor) const;
        Eigen::Matrix<Scalar, 6, 6> get_df_dx(const Vec3& x, const Vec3& theta) const;
};

#endif // RIGID_BODY_H_
