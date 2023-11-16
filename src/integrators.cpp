#include <Eigen/IterativeLinearSolvers>
#include <vector>
#include "integrators.hpp"
#include "linear_algebra.hpp"
#include "simulation.hpp"

void integrate_implicit_euler(const Simulation& simulation, PhysicsState* state, const EnergyAndDerivatives& f) {
    const Scalar h = simulation.TimeStep;
    const unsigned int nDoF = state->x.size();

    // Sparse Matrix creation
    // ----------------------------------------------------------------------------------
    SparseMat mass_matrix(nDoF, nDoF);
    const std::vector<Triplet> mass_matrix_triplets = compute_global_mass_matrix(simulation.simulables, *state);
    mass_matrix.setFromTriplets(mass_matrix_triplets.begin(), mass_matrix_triplets.end());
    SparseMat df_dx(nDoF, nDoF), df_dv(nDoF, nDoF);
    df_dx.setFromTriplets(f.df_dx_triplets.begin(), f.df_dx_triplets.end());
    // df_dv = h * df_dx;
    // df_dv.setFromTriplets(f.df_dv_triplets.begin(), f.df_dv_triplets.end());
    // ----------------------------------------------------------------------------------

    // Construct the system of equations
    // ----------------------------------------------------------------------------------
    Vec equation_vector = h * (f.force + h * df_dx * state->v);
    SparseMat equation_matrix = mass_matrix - h * df_dv - h * h * df_dx;
    handle_frozen_dof(simulation.frozen_dof, &equation_vector, &equation_matrix);
    // ----------------------------------------------------------------------------------

    // Solving the system of equations
    // ----------------------------------------------------------------------------------
    // Gradient conjugate solving method class
    Eigen::ConjugateGradient<Eigen::SparseMatrix<Scalar>> cg;
    cg.compute(equation_matrix);
    const Vec delta_v = cg.solve(equation_vector);
    // ----------------------------------------------------------------------------------

    // Update the state with the result
    // ----------------------------------------------------------------------------------
    state->v += delta_v;
    state->x += state->v * h;
}

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
