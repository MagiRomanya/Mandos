#ifndef HARD_CONSTRAINTS_H_
#define HARD_CONSTRAINTS_H_

#include "linear_algebra.hpp"
#include "physics_state.hpp"
#include "rigid_body.hpp"

struct RB_PointConstraint {
    RB_PointConstraint(unsigned int index, const RigidBody& rbA, const RigidBody& rbB, const Vec3& pA, const Vec3& pB) : index(index), pA(pA), pB(pB), rbA(rbA), rbB(rbB) {}
    static const unsigned int n_constraints = 3;
    const Vec3 pA, pB; // Joined points in rb frame of reference
    const RigidBody rbA, rbB;
    const unsigned int index;

    Vec3 evaluate_constraint(const PhysicsState& state) const;
    Eigen::Matrix<Scalar, 3, 12> evaluate_constraint_jacobian(const PhysicsState& state) const;

    void compute_constraint_and_jacobian(const PhysicsState & state, ConstraintsAndJacobians& c) const;
};

#endif // HARD_CONSTRAINTS_H_
