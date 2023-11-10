#ifndef RIGID_BODY_H_
#define RIGID_BODY_H_

#include "linear_algebra.hpp"
#include "physics_state.hpp"

struct RigidBodyParameters {
    Scalar k, L0;
};

struct RigidBody {
    RigidBody(unsigned int index, RigidBodyParameters param)
        : index(index), parameters(param) {}

    const unsigned int index;
    const RigidBodyParameters parameters;

    private:
        Scalar get_energy(const Vec3& x, const Vec3& theta);
        Eigen::Vector<Scalar, 6> get_force(const Vec3& x, const Vec3& theta);
        Eigen::Matrix<Scalar, 6, 6> get_df_dx(const Vec3& x, const Vec3& theta);
};

#endif // RIGID_BODY_H_
