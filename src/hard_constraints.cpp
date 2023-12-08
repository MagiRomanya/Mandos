#include "hard_constraints.hpp"
#include "utility_functions.hpp"

Vec3 RB_PointConstraint::evaluate_constraint(const PhysicsState& state) const {
    const Vec3 O1 = rbA.get_COM_position(state);
    const Vec3 O2 = rbB.get_COM_position(state);
    const Vec3 world_pos1 = rbA.compute_rotation_matrix(state.x) * pA + O1;
    const Vec3 world_pos2 = rbB.compute_rotation_matrix(state.x) * pB + O2;

    return world_pos1 - world_pos2; // = 0 when the constraint is met
}

Eigen::Matrix<Scalar, 3, 12> RB_PointConstraint::evaluate_constraint_jacobian(const PhysicsState& state) const {
    const Mat3 dc_dx = Mat3::Identity();
    const Mat3 dc_dtheta1 = - skew(rbA.compute_rotation_matrix(state.x) * pA);
    const Mat3 dc_dtheta2 = - skew(rbB.compute_rotation_matrix(state.x) * pB);
    Eigen::Matrix<Scalar, 3, 12> result;
    result << dc_dx, dc_dtheta1, dc_dx, dc_dtheta2;
    return result;
}

void RB_PointConstraint::compute_constraint_and_jacobian(const PhysicsState & state, ConstraintsAndJacobians& c) const {
    const unsigned int nDoF = state.get_nDoF();
    const unsigned int g_index = nDoF + index;

    c.constraints.segment<3>(index) = evaluate_constraint(state);

    const Eigen::Matrix<Scalar, 3, 12> jacobian = evaluate_constraint_jacobian(state);
    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 6; j++) {
            // Jacobian
            c.jacobian_triplets.emplace_back(g_index, rbA.index, jacobian(i,j));
            c.jacobian_triplets.emplace_back(g_index, rbB.index, jacobian(i,6+j));

            // Jacobian transposed
            c.jacobian_triplets.emplace_back(rbA.index, g_index, jacobian(j  ,i));
            c.jacobian_triplets.emplace_back(rbB.index, g_index, jacobian(6+j,i));
        }
    }
}
