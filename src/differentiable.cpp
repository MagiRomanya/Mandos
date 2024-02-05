#include <Eigen/IterativeLinearSolvers>
#include <cassert>
#include "differentiable.hpp"
#include "linear_algebra.hpp"
#include "utility_functions.hpp"

/**
 * The only energy gradients that depend on x0 are the ones which depend on v, which means inertial energies and some drag energies.
 */
void compute_partial_gradient_energy_partial_position(Scalar TimeStep, const Energies& energies, const PhysicsState& state, SparseMat& out) {
    std::vector<Triplet> dgradE_dx0_triplets;
    const Scalar one_over_h2 = 1.0f / (TimeStep * TimeStep);
    for (unsigned int i = 0; i < energies.linear_inertias.size(); i++) {
        LinearInertia e = energies.linear_inertias[i];
        const Mat3 dgradE_dx0 = - one_over_h2 * e.Mass;

        for (unsigned int i = 0; i < 3; i++) {
            for (unsigned int j = 0; j < 3; j++) {
                dgradE_dx0_triplets.emplace_back(e.p.index + i, e.p.index + j,  dgradE_dx0(i, j));
            }
        }
    }
    for (unsigned int i = 0; i < energies.rotational_inertias.size(); i++) {
        assert(false && "TODO: RIGID BODY NOT IMPLEMENTED");
    }
    // TODO: IGNORING DRAG!
    out.setFromTriplets(dgradE_dx0_triplets.begin(), dgradE_dx0_triplets.end());
}

/**
 * The only energy gradients that depend on v0 are the inertial energies and maybe some retarded friction model
 */
void compute_partial_gradient_energy_partial_velocity(Scalar TimeStep, const Energies& energies, const PhysicsState& state, SparseMat& out) {
    std::vector<Triplet> dgradE_dv0_triplets;
    const Scalar one_over_h2 = 1.0f / (TimeStep * TimeStep);
    for (unsigned int i = 0; i < energies.linear_inertias.size(); i++) {
        LinearInertia e = energies.linear_inertias[i];
        const Mat3 dgradE_dv0 = - TimeStep * one_over_h2 * e.Mass;

        for (unsigned int i = 0; i < 3; i++) {
            for (unsigned int j = 0; j < 3; j++) {
                dgradE_dv0_triplets.emplace_back(e.p.index + i, e.p.index + j,  dgradE_dv0(i, j));
            }
        }
    }
    for (unsigned int i = 0; i < energies.rotational_inertias.size(); i++) {
        assert(false && "TODO: RIGID BODY NOT IMPLEMENTED");
    }

    // TODO: IGNORING DRAG!
    out.setFromTriplets(dgradE_dv0_triplets.begin(), dgradE_dv0_triplets.end());
}

Vec compute_loss_function_gradient_backpropagation(const Simulation& simulation,
                                                   const std::vector<PhysicsState>& trajectory,
                                                   const LossFunctionAndDerivatives& loss,
                                                   const Mat& dx0_dp, const Mat& dv0_dp)
{
    const unsigned int nParameters = loss.loss_parameter_partial_derivative.size();
    const unsigned int nDoF = trajectory.at(0).x.size();
    const unsigned int nStates = trajectory.size();

    // Initialize the loss function gradients dg_dp, dg_dx and dg_dv
    // -------------------------------------------------------------------------
    Vec loss_gradient = loss.loss_parameter_partial_derivative;
    Vec loss_position_gradient = loss.loss_position_partial_derivative[nStates-1];
    Vec loss_velocity_gradient = loss.loss_velocity_partial_derivative[nStates-1];

    // Backward loop
    // -------------------------------------------------------------------------
    const Scalar one_over_h = 1.0f / simulation.TimeStep;
    for (int i = nStates-2; i >= 0; i--) {
        const PhysicsState& state = trajectory[i+1];

        // Compute the hessian matrix and other sparse matrices
        // -------------------------------------------------------------------------
        EnergyAndDerivatives f(nDoF);
        compute_energy_and_derivatives(simulation.TimeStep, simulation.energies, state, state, f);
        SparseMat hessian(nDoF,nDoF), dgradE_dx0(nDoF,nDoF), dgradE_dv0(nDoF,nDoF);
        Mat dgradE_dp = Mat::Zero(nDoF, nParameters);
        hessian.setFromTriplets(f.hessian_triplets.begin(), f.hessian_triplets.end());
        compute_partial_gradient_energy_partial_position(simulation.TimeStep, simulation.energies, state, dgradE_dx0);
        compute_partial_gradient_energy_partial_velocity(simulation.TimeStep, simulation.energies, state, dgradE_dv0);

        // Solve the linear system
        // -------------------------------------------------------------------------
        Eigen::ConjugateGradient<SparseMat> cg;
        cg.compute(hessian);
        const Vec equation_vector = - one_over_h * loss_velocity_gradient - loss_position_gradient;
        Vec alpha = cg.solve(equation_vector);

        // Update the loss function gradients
        // -------------------------------------------------------------------------
        loss_position_gradient = loss.loss_position_partial_derivative[i] - one_over_h * loss_velocity_gradient + dgradE_dx0 * alpha;
        loss_velocity_gradient = loss.loss_velocity_partial_derivative[i] + dgradE_dv0 * alpha;
        loss_gradient += dgradE_dp * alpha;
    }

    // Add the initial conditions term
    // -------------------------------------------------------------------------
    loss_gradient += loss_position_gradient * dx0_dp + loss_velocity_gradient * dv0_dp;

    return loss_gradient;
}
