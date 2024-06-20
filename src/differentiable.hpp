#ifndef DIFFERENTIABLE_H_
#define DIFFERENTIABLE_H_

#include "simulation.hpp"

struct LossFunctionAndDerivatives {
    Scalar loss;                                       // g
    std::vector<Vec> loss_position_partial_derivative; // dg_dx
    std::vector<Vec> loss_velocity_partial_derivative; // dg_dv
    Vec loss_parameter_partial_derivative;             // dg_dp
};

Vec compute_loss_function_gradient_backpropagation(const Simulation& simulation,
                                                   const std::vector<PhysicsState>& trajectory,
                                                   const LossFunctionAndDerivatives& loss,
                                                   const Mat& dx0_dp, const Mat& dv0_dp,
                                                   const unsigned int maxIterations = 0);

Vec compute_loss_function_gradient_backpropagation_control(const Simulation& simulation,
                                                           const std::vector<PhysicsState>& trajectory,
                                                           const LossFunctionAndDerivatives& loss,
                                                           const Mat& dx0_dp, const Mat& dv0_dp,
                                                           const unsigned int maxIterations);

Vec compute_loss_function_gradient_backpropagation_1_step_velocity(const Simulation& simulation,
                                                                   const PhysicsState state0,
                                                                   const PhysicsState state1,
                                                                   const Vec dg_dphi,
                                                                   const Vec dg_dphi_dot);

Vec compute_loss_function_gradient_backpropagation_1_step_position(const Simulation& simulation,
                                                                   const PhysicsState state0,
                                                                   const PhysicsState state1,
                                                                   const Vec dg_dphi,
                                                                   const Vec dg_dphi_dot);

#endif // DIFFERENTIABLE_H_
