#ifndef MANDOS_INTEGRATORS_H_
#define MANDOS_INTEGRATORS_H_

#include <Mandos/linear_algebra.hpp>
#include <Mandos/simulation.hpp>

namespace mandos
{

void handle_frozen_dof(const std::vector<unsigned int>& frozen_dof, Vec* eq_vec, SparseMat* eq_mat);

void integrate_implicit_euler(const Simulation& simulation,
                              const PhysicsState& state,
                              const EnergyAndDerivatives& f,
                              Vec& dx);

}  // namespace mandos

#endif  // MANDOS_INTEGRATORS_H_
