#ifndef INTEGRATORS_H_
#define INTEGRATORS_H_

#include "linear_algebra.hpp"
#include "simulation.hpp"

void handle_frozen_dof(const std::vector<unsigned int>& frozen_dof, Vec* eq_vec, SparseMat* eq_mat);

void integrate_implicit_euler(const Simulation& simulation, const PhysicsState& state, const EnergyAndDerivatives& f, const ConstraintsAndJacobians& c, Vec& dx);

#endif // INTEGRATORS_H_
