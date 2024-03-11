#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#include <Eigen/Dense> // For inverse
#include "differentiable.hpp"
#include "inertia_energies.hpp"
#include "linear_algebra.hpp"
#include "physics_state.hpp"
#include "rigid_body.hpp"
#include "simulation.hpp"
#include "utility_functions.hpp"

static const Scalar threshold2 = 1e-5;

inline Mat3 compute_global_axis_angle_jacobian(const Vec3& phi) {
    const Scalar angle = phi.norm();
    if (angle < threshold2) return Mat3::Identity();
    const Vec3 axis = phi / angle;
    const Scalar half_angle = 0.5 * angle;
    const Mat3 axisaxisT = axis * axis.transpose();
    return half_angle / tan(half_angle) * (Mat3::Identity() - axisaxisT)
        + axisaxisT
        - skew(0.5 * phi);
}

inline void compute_axis_angle_jacobian_parts(const Vec3& phi, Mat3& A, Mat3& B) {
    const Mat3 phi_phiT = phi * phi.transpose();
    A = compute_rotation_matrix_rodrigues(phi) - Mat3::Identity() + phi_phiT;
    B = skew(phi) + phi_phiT;
}

inline Mat3 compute_global_to_local_axis_angle_jacobian(const Vec3& phi) {
    if (phi.squaredNorm() < threshold2) return Mat3::Identity();

    Mat3 A, B;
    compute_axis_angle_jacobian_parts(phi, A, B);
    return A.inverse() * B;
}

inline Mat3 compute_local_to_global_axis_angle_jacobian(const Vec3& phi) {
    if (phi.squaredNorm() < threshold2) return Mat3::Identity();

    Mat3 A, B;
    compute_axis_angle_jacobian_parts(phi, A, B);
    return B.inverse() * A;
}


/**
 * Compute the derivative of a rotation matrix with respect to local axis angle, evaluated at theta.
 */
Eigen::Matrix<Scalar,3,9> dvecR_dtheta_local(const Vec3& theta) {
    const Mat3 R = compute_rotation_matrix_rodrigues(theta);
    return vectorized_levi_civita() * block_matrix(R);
}

/**
 * Compute the derivative of a rotation matrix with respect to the global axis angle theta.
 */
inline Eigen::Matrix<Scalar,3,9> dvecR_dtheta_global(const Vec3& theta) {
    const Mat3 jac = compute_local_to_global_axis_angle_jacobian(theta);
    const Mat3 R = compute_rotation_matrix_rodrigues(theta);
    return jac.transpose() * vectorized_levi_civita() * block_matrix(R);
}

Mat3 rotation_inertia_energy_hessian(const Mat3& J_inertia_tensor, const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    const Mat3 R = compute_rotation_matrix_rodrigues(theta);
    const Mat3 R0 = compute_rotation_matrix_rodrigues(theta0);
    const Mat3 R0old = compute_rotation_matrix_rodrigues(theta0 - omega0 * TimeStep);
    const Mat3 R_guess = (R0 + (R0 - R0old)); // x0 + h* v0

    const Mat3 rot_inertia = R * J_inertia_tensor * R_guess.transpose();
    const Mat3 S = (rot_inertia + rot_inertia.transpose()) / 2; // Exact hessian
    // const Mat3 S = R * J_inertia_tensor * R.transpose(); // Linear approximation
    const Scalar h2 = TimeStep * TimeStep;

    const Mat3 hessian = 1.0 / h2 * (S.trace() * Mat3::Identity() - S);
    return hessian;
}

inline Mat3 rotation_inertia_dgradE_dtheta(const Mat3& J_inertia_tensor, const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    const Mat3 R0 = compute_rotation_matrix_rodrigues(theta0);
    const Mat3 R0old = compute_rotation_matrix_rodrigues(theta0 - TimeStep * omega0);
    const Mat3 Rguess = (R0 + (R0 - R0old)); // x0 + h* v0

    const Scalar h2 = TimeStep * TimeStep;
    const Eigen::Matrix<Scalar,3,9> vLeviCivita = vectorized_levi_civita();
    const Eigen::Matrix<Scalar,3,9> dvecR_dtheta = dvecR_dtheta_global(theta);

    const Eigen::Matrix<Scalar,9,3> dvecRMR_guess_dtheta = block_matrix<3,3>(Rguess * J_inertia_tensor) * dvecR_dtheta.transpose();;
    const Eigen::Matrix<Scalar,9,3> dvecAdtheta = 0.5 * (dvecRMR_guess_dtheta - transpose_vectorized_matrix(dvecRMR_guess_dtheta));

    Mat3 H = 1.0 / h2 * vLeviCivita * dvecAdtheta;
    return H;
}

inline Mat3 rotation_inertia_dgradE_dtheta0(const Mat3& J_inertia_tensor, const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    const Mat3 R = compute_rotation_matrix_rodrigues(theta);

    const Scalar h2 = TimeStep * TimeStep;
    const Eigen::Matrix<Scalar,3,9> vLeviCivita = vectorized_levi_civita();
    const Eigen::Matrix<Scalar,3,9> dvecRguess_dtheta0 = 2 * dvecR_dtheta_global(theta0) - dvecR_dtheta_global(theta0 - omega0 * TimeStep);

    const Eigen::Matrix<Scalar,9,3> dvecRMR_guess_dtheta0 = block_matrix<3,3>(R * J_inertia_tensor) * dvecRguess_dtheta0.transpose();;
    const Eigen::Matrix<Scalar,9,3> dvecAdtheta0 = 0.5 * (transpose_vectorized_matrix(dvecRMR_guess_dtheta0) - dvecRMR_guess_dtheta0);

    Mat3 H = 1.0 / h2 * vLeviCivita * dvecAdtheta0;
    return H;
}

inline Mat3 rotation_inertia_dgradE_domega0(const Mat3& J_inertia_tensor, const Vec3& theta, const Vec3& theta0, const Vec3& omega0, Scalar TimeStep) {
    const Mat3 R = compute_rotation_matrix_rodrigues(theta);

    const Scalar h2 = TimeStep * TimeStep;
    const Eigen::Matrix<Scalar,3,9> vLeviCivita = vectorized_levi_civita();
    const Eigen::Matrix<Scalar,3,9> dvecRguess_domega0 = TimeStep * dvecR_dtheta_global(theta0 - omega0 * TimeStep);

    const Eigen::Matrix<Scalar,9,3> dvecRMR_guess_domega0 = block_matrix<3,3>(R * J_inertia_tensor) * dvecRguess_domega0.transpose();;
    const Eigen::Matrix<Scalar,9,3> dvecAdtheta0 = 0.5 * (transpose_vectorized_matrix(dvecRMR_guess_domega0) - dvecRMR_guess_domega0);

    Mat3 H = 1.0 / h2 * vLeviCivita * dvecAdtheta0;
    return H;
}

// #define DIFF_LOCAL_RB
void compute_gradient_partial_derivatives(Scalar TimeStep, const Energies& energies, const PhysicsState& state, const PhysicsState& state0, SparseMat& dgradE_dx, SparseMat& dgradE_dx0, SparseMat& dgradE_dv0) {
    const unsigned int nDoF = state.get_nDoF();
    EnergyAndDerivatives f(nDoF);
    std::vector<Triplet> dgradE_dx0_triplets;
    std::vector<Triplet> dgradE_dv0_triplets;

    // Linear inertias
    // ----------------------------------------------------------------------------------------------------
    const Scalar one_over_h2 = 1.0f / (TimeStep * TimeStep);
    for (unsigned int i = 0; i < energies.linear_inertias.size(); i++) {
        LinearInertia e = energies.linear_inertias[i];
        const Mat3 dgradE_dx = one_over_h2 * e.Mass;
        const Mat3 dgradE_dx0 = - dgradE_dx;
        const Mat3 dgradE_dv0 = - TimeStep * dgradE_dx;

        for (unsigned int i = 0; i < 3; i++) {
            for (unsigned int j = 0; j < 3; j++) {
                f.hessian_triplets.emplace_back(e.p.index + i, e.p.index + j,  dgradE_dx(i, j));
                dgradE_dx0_triplets.emplace_back(e.p.index + i, e.p.index + j,  dgradE_dx0(i, j));
                dgradE_dv0_triplets.emplace_back(e.p.index + i, e.p.index + j,  dgradE_dv0(i, j));
            }
        }
    }

    // Rotational inertias
    // ----------------------------------------------------------------------------------------------------
    for (unsigned int i = 0; i < energies.rotational_inertias.size(); i++) {
        const RotationalInertia& e = energies.rotational_inertias[i];
        const Mat3 J_inertia_tensor = e.rb.J_inertia_tensor0;
        const Vec3 theta = e.rb.get_axis_angle(state.x);
        const Vec3 theta0 = e.rb.get_axis_angle(state0.x);
        const Vec3 omega0 = e.rb.get_axis_angle(state0.v);
#ifdef DIFF_LOCAL_RB
        const Mat3 dgradE_dtheta = rotation_inertia_energy_hessian(J_inertia_tensor, theta, theta0, omega0, TimeStep);
#else
        const Mat3 dgradE_dtheta = rotation_inertia_dgradE_dtheta(J_inertia_tensor, theta, theta0, omega0, TimeStep);
#endif
        const Mat3 dgradE_dtheta0 = rotation_inertia_dgradE_dtheta0(J_inertia_tensor, theta, theta0, omega0, TimeStep);
        const Mat3 dgradE_domega0 = rotation_inertia_dgradE_domega0(J_inertia_tensor, theta, theta0, omega0, TimeStep);

        for (unsigned int i = 0; i < 3; i++) {
            for (unsigned int j = 0; j < 3; j++) {
                f.hessian_triplets.emplace_back(e.rb.index + 3 + i, e.rb.index + 3 + j,  dgradE_dtheta(i, j));
                dgradE_dx0_triplets.emplace_back(e.rb.index + 3 + i, e.rb.index + 3 + j,  dgradE_dtheta0(i, j));
                dgradE_dv0_triplets.emplace_back(e.rb.index + 3 + i, e.rb.index + 3 + j,  dgradE_domega0(i, j));
            }
        }
    }


    // Potential energies
    // ---------------------------------------------------------------------------------
#define MAT(type, name) X(std::vector<FEM_Element3D<type>>, fem_elements_##name)
#define X(type, energy) \
    for (size_t i = 0; i < energies.energy.size(); i++) { \
        energies.energy[i].compute_energy_and_derivatives(TimeStep, state, f); \
    }
    POTENTIAL_ENERGY_MEMBERS
#undef X
#undef MAT


    // NOTE: IGNORING POTENTIALS WITH VELOCITY DEPENDENCE

    dgradE_dx.setFromTriplets(f.hessian_triplets.begin(), f.hessian_triplets.end());
    dgradE_dx0.setFromTriplets(dgradE_dx0_triplets.begin(), dgradE_dx0_triplets.end());
    dgradE_dv0.setFromTriplets(dgradE_dv0_triplets.begin(), dgradE_dv0_triplets.end());
}

void compute_simulation_global_to_local_jacobian(const Simulables& simulables, const PhysicsState state, SparseMat& ptm) {
    std::vector<Triplet> jac_triplets;
    for (unsigned int i = 0; i < simulables.particles.size(); i++) {
        const unsigned int idx = simulables.particles[i].index;
        for (unsigned int j = 0; j < 3; j++) jac_triplets.emplace_back(idx + j, idx + j, 1); // Diagonal matrix
    }
    for (unsigned int i = 0; i < simulables.rigid_bodies.size(); i++) {
        const RigidBody& rb = simulables.rigid_bodies[i];
        const Vec3 theta = rb.get_axis_angle(state.x);
        const Mat3 jac = compute_global_to_local_axis_angle_jacobian(theta);
        // const Mat3 jac = compute_global_axis_angle_jacobian(theta);

        const unsigned int idx = rb.index + 3;

        for (unsigned int j = 0; j < 3; j++) {
            jac_triplets.emplace_back(rb.index + j, rb.index + j, 1); // Translation diagonal
            for (unsigned int k = 0; k < 3; k++) {
                jac_triplets.emplace_back(idx + j, idx + k, jac(j, k));
            }
        }
    }
    ptm.setFromTriplets(jac_triplets.begin(), jac_triplets.end());
}


Vec compute_loss_function_gradient_backpropagation(const Simulation& simulation,
                                                   const std::vector<PhysicsState>& trajectory,
                                                   const LossFunctionAndDerivatives& loss,
                                                   const Mat& dx0_dp, const Mat& dv0_dp)
{
    const unsigned int nParameters = loss.loss_parameter_partial_derivative.size();
    const unsigned int nDoF = trajectory.at(0).x.size();
    const unsigned int nStates = trajectory.size();
    const unsigned int nSteps = nStates - 1;

    // Initialize the loss function gradients dg_dp, dg_dx and dg_dv
    // -------------------------------------------------------------------------
    Vec loss_gradient = loss.loss_parameter_partial_derivative;
    Vec loss_position_gradient = loss.loss_position_partial_derivative[nSteps];
    Vec loss_velocity_gradient = loss.loss_velocity_partial_derivative[nSteps];

    // Backward loop
    // -------------------------------------------------------------------------
    const Scalar one_over_h = 1.0f / simulation.TimeStep;
    for (int i = nSteps - 1; i >= 0; i--) {
        const PhysicsState& state = trajectory[i+1];
        const PhysicsState& state0 = trajectory[i];

        // Compute the hessian matrix and other sparse matrices
        // -------------------------------------------------------------------------

        SparseMat hessian(nDoF,nDoF), dgradE_dx0(nDoF,nDoF), dgradE_dv0(nDoF,nDoF);
        compute_gradient_partial_derivatives(simulation.TimeStep, simulation.energies, state, state0, hessian, dgradE_dx0, dgradE_dv0);

        // NOTE: This is not always zero
        Mat dgradE_dp = Mat::Zero(nDoF, nParameters);

        // Solve the linear system
        // -------------------------------------------------------------------------
#ifdef DIFF_LOCAL_RB
        Eigen::ConjugateGradient<SparseMat> solver;
        SparseMat jac_global_to_local(nDoF, nDoF);
        compute_simulation_global_to_local_jacobian(simulation.simulables, state, jac_global_to_local);
        const Vec equation_vector = - (loss_position_gradient + one_over_h * loss_velocity_gradient).transpose() * jac_global_to_local;
#else
        const Vec equation_vector = - (loss_position_gradient + one_over_h * loss_velocity_gradient);
        Eigen::SparseLU<SparseMat> solver;
#endif
        solver.compute(hessian.transpose());
        Vec alpha = solver.solve(equation_vector);

        // Update the loss function gradients
        // -------------------------------------------------------------------------
        loss_position_gradient = (loss.loss_position_partial_derivative[i].transpose() - one_over_h * loss_velocity_gradient.transpose() + alpha.transpose() * dgradE_dx0);
        loss_velocity_gradient = (loss.loss_velocity_partial_derivative[i].transpose() + alpha.transpose() * dgradE_dv0);

        loss_gradient += alpha.transpose() * dgradE_dp;
    }

    // Add the initial conditions term
    // -------------------------------------------------------------------------
    loss_gradient += loss_position_gradient.transpose() * dx0_dp + loss_velocity_gradient.transpose() * dv0_dp;

    return loss_gradient;
}

Vec compute_loss_function_gradient_backpropagation_1_step_velocity(const Simulation& simulation,
                                                                   const PhysicsState state0,
                                                                   const PhysicsState state1,
                                                                   const Vec dg_dphi,
                                                                   const Vec dg_dphi_dot)
{
    const unsigned int nDoF = state0.x.size();

    SparseMat hessian(nDoF,nDoF), dgradE_dx0(nDoF,nDoF), dgradE_dv0(nDoF,nDoF);
    compute_gradient_partial_derivatives(simulation.TimeStep, simulation.energies, state1, state0, hessian, dgradE_dx0, dgradE_dv0);

    const Mat dgradE_dphi = hessian.toDense();
    const Mat dgradE_dphi0 = dgradE_dx0.toDense();
    const Mat dgradE_dphi_dot0 = dgradE_dv0.toDense();

    SparseMat ptm(nDoF, nDoF);
    compute_simulation_global_to_local_jacobian(simulation.simulables, state1, ptm);

    Eigen::ConjugateGradient<Mat> cg;
    cg.compute(-dgradE_dphi);
    const Vec eq_vec = dg_dphi.transpose() + 1.0 / simulation.TimeStep * dg_dphi_dot.transpose() * ptm;
    // const Vec eq_vec = dg_dphi.transpose() + 1.0 / simulation.TimeStep * dg_dphi_dot.transpose();
    Vec adjoint = cg.solve(eq_vec);

    const Vec dgdp_adj = adjoint.transpose() * dgradE_dphi_dot0;

    return dgdp_adj;
}

Vec compute_loss_function_gradient_backpropagation_1_step_position(const Simulation& simulation,
                                                                   const PhysicsState state0,
                                                                   const PhysicsState state1,
                                                                   const Vec dg_dphi,
                                                                   const Vec dg_dphi_dot)
{
    const unsigned int nDoF = state0.x.size();

    SparseMat hessian(nDoF,nDoF), dgradE_dx0(nDoF,nDoF), dgradE_dv0(nDoF,nDoF);
    compute_gradient_partial_derivatives(simulation.TimeStep, simulation.energies, state1, state0, hessian, dgradE_dx0, dgradE_dv0);

    const Mat dgradE_dphi = hessian.toDense();
    const Mat dgradE_dphi0 = dgradE_dx0.toDense();
    const Mat dgradE_dphi_dot0 = dgradE_dv0.toDense();

    SparseMat ptm(nDoF, nDoF);
    compute_simulation_global_to_local_jacobian(simulation.simulables, state1, ptm);

    Eigen::ConjugateGradient<Mat> cg;
    cg.compute(-dgradE_dphi);
    const Vec eq_vec = dg_dphi.transpose()*ptm + 1.0 / simulation.TimeStep * dg_dphi_dot.transpose() * ptm;
    // const Vec eq_vec = dg_dphi;
    Vec adjoint = cg.solve(eq_vec);

    const Mat dphi_dphi0 = - dgradE_dphi.inverse() * dgradE_dphi0;

    const Vec dgdp_adj = - 1.0 / simulation.TimeStep * dg_dphi_dot.transpose() + adjoint.transpose() * dgradE_dphi0;

    const Vec dgdp = dg_dphi.transpose() * dphi_dphi0;

    // return dgdp;
    return dgdp_adj;
}
