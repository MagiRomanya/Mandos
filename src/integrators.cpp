#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include <iostream>
#include <vector>

#include "integrators.hpp"
#include "linear_algebra.hpp"
#include "particle_rigid_body_copuling.hpp"
#include "physics_state.hpp"
#include "simulation.hpp"

void integrate_implicit_euler(const Simulation& simulation, const PhysicsState& state, const EnergyAndDerivatives& f, Vec& dx) {
    const unsigned int nDoF = state.get_nDoF();

    SparseMat equation_matrix(nDoF, nDoF);
    Vec equation_vector = - f.gradient;
    equation_matrix.setFromTriplets(f.hessian_triplets.begin(), f.hessian_triplets.end());

    handle_frozen_dof(simulation.frozen_dof, &equation_vector, &equation_matrix);

    SparseMat copuling_jacobian;
    compute_copuling_jacobian(simulation.copulings, state, copuling_jacobian);
    std::cout << copuling_jacobian.toDense() << std::endl;
    const SparseMat copuling_jacobian_t = copuling_jacobian.transpose();

    equation_matrix = copuling_jacobian_t * equation_matrix * copuling_jacobian;
    equation_vector = copuling_jacobian_t * equation_vector;
    // ----------------------------------------------------------------------------------

    // Solving the system of equations
    // ----------------------------------------------------------------------------------
    // Gradient conjugate solving method class
    Eigen::ConjugateGradient<Eigen::SparseMatrix<Scalar>> cg;
    cg.compute(equation_matrix);
    dx = copuling_jacobian * cg.solve(equation_vector);
}

#ifdef ENABLE_LAGRANGE_MULTIPLIER_CONSTRAINTS
void integrate_implicit_euler(const Simulation& simulation, const PhysicsState& state, const EnergyAndDerivatives& f, const ConstraintsAndJacobians& c, Vec& dx) {
    const unsigned int nDoF = state.get_nDoF();
    const unsigned int n_constraints = c.get_n_constraints();

    // ----------------------------------------------------------------------------------
    std::vector<Triplet> equation_matrix_triplets;
    equation_matrix_triplets.insert(equation_matrix_triplets.end(), f.hessian_triplets.begin(), f.hessian_triplets.end());
    equation_matrix_triplets.insert(equation_matrix_triplets.end(), c.jacobian_triplets.begin(), c.jacobian_triplets.end());
    // ----------------------------------------------------------------------------------

    // Construct the system of equations
    // / H  J^T \ /  dx  \ = /-GradE\
    // \ J   0  / \lambda/ = \-C(x0)/
    // ----------------------------------------------------------------------------------
    Vec equation_vector(nDoF + n_constraints);
    equation_vector << - f.gradient, -c.constraints;
    SparseMat equation_matrix(nDoF + n_constraints, nDoF + n_constraints);
    equation_matrix.setFromTriplets(equation_matrix_triplets.begin(), equation_matrix_triplets.end());

    handle_frozen_dof(simulation.frozen_dof, &equation_vector, &equation_matrix);
    // ----------------------------------------------------------------------------------

    // Solving the system of equations
    // ----------------------------------------------------------------------------------
    // Gradient conjugate solving method class
    Eigen::ConjugateGradient<Eigen::SparseMatrix<Scalar>> cg;
    cg.compute(equation_matrix);
    dx = cg.solve(equation_vector);
}
#endif // ENABLE_LAGRANGE_MULTIPLIER_CONSTRAINTS

struct FrozenDoFPredicate {
    FrozenDoFPredicate(const std::vector<unsigned int>& frozen_dof) : frozen_dof(frozen_dof) {}

    bool operator() (const Eigen::Index& row, const Eigen::Index& col, const double& value) const {
        // Keep elements in the diagonal and outside the dof column and row
        for (unsigned int i = 0; i < frozen_dof.size(); i++) {
            unsigned int dof = frozen_dof[i];
            if (row==col) return true;
            else if ((row==dof) or (col==dof)) return false;
        }
        return true;
    }
    private:
        const std::vector<unsigned int>& frozen_dof;
};

void handle_frozen_dof(const std::vector<unsigned int>& frozen_dof, Vec* eq_vec, SparseMat* eq_mat) {
    // Eliminate non zeros from the rows and columns
    (*eq_mat).prune(FrozenDoFPredicate(frozen_dof));
    for (unsigned int i = 0; i < frozen_dof.size(); i++) {
        (*eq_vec)[frozen_dof[i]] = 0.0;
    }
}

void handle_frozen_dof(const std::vector<unsigned int>& frozen_dof, SparseMat* mat) {
    // Eliminate non zeros from the rows and columns
    (*mat).prune(FrozenDoFPredicate(frozen_dof));
}

void handle_frozen_dof(const std::vector<unsigned int>& frozen_dof, Vec* vec) {
    for (unsigned int i = 0; i < frozen_dof.size(); i++) {
        (*vec)[frozen_dof[i]] = 0.0;
    }
}
