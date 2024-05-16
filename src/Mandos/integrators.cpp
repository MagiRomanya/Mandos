#include <spdlog/spdlog.h>
#include <thread>
#include <vector>

#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

#include <Mandos/linear_algebra.hpp>

#include <Mandos/integrators.hpp>
#include <Mandos/particle_rigid_body_coupling.hpp>
#include <Mandos/physics_state.hpp>
#include <Mandos/simulation.hpp>
#include <Mandos/utility_functions.hpp>

namespace mandos
{

void integrate_implicit_euler(const Simulation& simulation,
                              const PhysicsState& state,
                              const EnergyAndDerivatives& f,
                              Vec& dx)
{
    const unsigned int nDoF = state.get_nDoF();

    SparseMat equation_matrix(nDoF, nDoF);
    Vec equation_vector = -f.gradient;
    equation_matrix.setFromTriplets(f.hessian_triplets.begin(), f.hessian_triplets.end());

    handle_frozen_dof(simulation.frozen_dof, &equation_vector, &equation_matrix);

    SparseMat coupling_jacobian;
    compute_coupling_jacobian(simulation.couplings, state, coupling_jacobian);
    const SparseMat coupling_jacobian_t = coupling_jacobian.transpose();

    equation_matrix = coupling_jacobian_t * equation_matrix * coupling_jacobian;
    equation_vector = coupling_jacobian_t * equation_vector;
    if (equation_vector.hasNaN()) {
        spdlog::warn("WARNING::INTEGRATE_IMPLICIT_EULER: Equation vector has NaN");
    }
    // ----------------------------------------------------------------------------------

    // Solving the system of equations
    // ----------------------------------------------------------------------------------
    Eigen::setNbThreads(std::thread::hardware_concurrency());
    Eigen::ConjugateGradient<SparseMat, Eigen::Lower | Eigen::Upper> solver;
    // Eigen::SparseLU<SparseMat> solver;
    const Scalar tol = 1e-9;
    double integration_solve_time = 0;
    {
        // TODO Use tracy
        // Clock clock(integration_solve_time);
        // solver.setTolerance(tol);
        solver.compute(equation_matrix);
        dx = coupling_jacobian * solver.solve(equation_vector);
    }

    if (dx.hasNaN()) {
        spdlog::info("WARNING::INTEGRATE_IMPLICIT_EULER: The dx step has NaN");
    }
}

struct FrozenDoFPredicate {
    FrozenDoFPredicate(const std::vector<unsigned int>& frozen_dof)
        : frozen_dof(frozen_dof)
    {
    }

    bool operator()(const Eigen::Index& row, const Eigen::Index& col, const double& value) const
    {
        // Keep elements in the diagonal and outside the dof column and row
        for (unsigned int i = 0; i < frozen_dof.size(); i++) {
            unsigned int dof = frozen_dof[i];
            if (row == col)
                return true;
            else if ((row == dof) or (col == dof))
                return false;
        }
        return true;
    }

private:
    const std::vector<unsigned int>& frozen_dof;
};

void handle_frozen_dof(const std::vector<unsigned int>& frozen_dof, Vec* eq_vec, SparseMat* eq_mat)
{
    // Eliminate non zeros from the rows and columns
    (*eq_mat).prune(FrozenDoFPredicate(frozen_dof));
    for (unsigned int i = 0; i < frozen_dof.size(); i++) {
        (*eq_vec)[frozen_dof[i]] = 0.0;
    }
}

void handle_frozen_dof(const std::vector<unsigned int>& frozen_dof, SparseMat* mat)
{
    // Eliminate non zeros from the rows and columns
    (*mat).prune(FrozenDoFPredicate(frozen_dof));
}

void handle_frozen_dof(const std::vector<unsigned int>& frozen_dof, Vec* vec)
{
    for (unsigned int i = 0; i < frozen_dof.size(); i++) {
        (*vec)[frozen_dof[i]] = 0.0;
    }
}

}  // namespace mandos
