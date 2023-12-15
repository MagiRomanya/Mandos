#ifndef HARD_CONSTRAINTS_H_
#define HARD_CONSTRAINTS_H_

// Toggle this define to enable/disable hard constraints using lagrange multipliers
// #define ENABLE_LAGRANGE_MULTIPLIER_CONSTRAINTS

#ifdef ENABLE_LAGRANGE_MULTIPLIER_CONSTRAINTS
#include "linear_algebra.hpp"
#include "physics_state.hpp"
#include "rigid_body.hpp"

struct ConstraintsAndJacobians {
    ConstraintsAndJacobians(unsigned int nConstraints) {
        constraints.setZero(nConstraints);
        jacobian_triplets.clear();
    }
    Vec constraints;
    std::vector<Triplet> jacobian_triplets;

    inline void increase_n_constraints(unsigned int increment_constraints) {
        constraints.conservativeResize(constraints.size() + increment_constraints);
    }

    inline unsigned int get_n_constraints() const { return constraints.size(); }
};

struct RB_PointConstraint {
    RB_PointConstraint(unsigned int index, const RigidBody& rbA, const RigidBody& rbB, const Vec3& pA, const Vec3& pB) : index(index), pA(pA), pB(pB), rbA(rbA), rbB(rbB) {}
    const Vec3 pA, pB; // Joined points in rb frame of reference
    const RigidBody rbA, rbB;
    const unsigned int index;

    Vec3 evaluate_constraint(const PhysicsState& state) const;
    Eigen::Matrix<Scalar, 3, 12> evaluate_constraint_jacobian(const PhysicsState& state) const;

    void compute_constraint_and_jacobian(const PhysicsState & state, ConstraintsAndJacobians& c) const;
};

struct HardConstraints {
    std::vector<RB_PointConstraint> rb_point_constraints;
};

void compute_constraints_and_jacobians(const HardConstraints& c, const PhysicsState& state, ConstraintsAndJacobians& out);

inline unsigned int count_constraints(const HardConstraints& constraints) {
    unsigned int nConstraints = 0;
    nConstraints += constraints.rb_point_constraints.size() * 3;
    return nConstraints;
}
#endif // ENABLE_LAGRANGE_MULTIPLIER_CONSTRAINTS

#endif // HARD_CONSTRAINTS_H_
