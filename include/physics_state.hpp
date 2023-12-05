#ifndef PHYSICSSTATE_H_
#define PHYSICSSTATE_H_

#include <vector>
#include <Eigen/Geometry>

#include "linear_algebra.hpp"

struct PhysicsState {
    Vec x;
    Vec x_old;

    void add_size(unsigned int increment_dof) {
        x.conservativeResize(x.size() + increment_dof);
        x_old.conservativeResize(x_old.size() + increment_dof);
    }
};

struct EnergyAndDerivatives {
    EnergyAndDerivatives(unsigned int nDoF) {
        energy = 0;
        gradient.setZero(nDoF);
    }
    // Container
    Scalar energy;
    Vec gradient;
    std::vector<Triplet> hessian_triplets;
};

inline void set_linear_velocity(PhysicsState& state, Scalar TimeStep, unsigned int index, const Vec3& v) {
    state.x_old.segment<3>(index) = state.x.segment<3>(index) - TimeStep * v;
}

inline void set_angular_velocity(PhysicsState& state, Scalar TimeStep, unsigned int index, const Vec3& omega) {
    const Vec3 theta = state.x.segment<3>(index);
    const Scalar angle = theta.norm();

    if (theta.norm() == 0) {
        state.x_old.segment<3>(index) = -omega;
        return;
    }
    const Vec3 axis = theta / theta.norm();
    typedef Eigen::Quaternion<Scalar> Quat;
    const Quat q = Quat(Eigen::AngleAxis<Scalar>(angle, axis));
    const Quat q_omega = Quat(0, -0.5 * omega);
    const Quat q_dot = q_omega * q;
    const Quat q_new = Quat(q.w() + q_dot.w(), q.vec() + q_dot.vec());
    const Eigen::AngleAxis<Scalar> angle_axis(q_new);

    state.x_old.segment<3>(index) = angle_axis.angle() * angle_axis.axis();
}


#endif // PHYSICSSTATE_H_
