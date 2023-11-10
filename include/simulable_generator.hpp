#ifndef SIMULABLE_GENERATOR_H_
#define SIMULABLE_GENERATOR_H_

#include "simulation.hpp"
#include "linear_algebra.hpp"
#include "physics_state.hpp"

void generate_mass_spring(Simulation& simulation,
                          const std::vector<Scalar>& vertices,
                          const std::vector<unsigned int>& indices,
                          Scalar node_mass,
                          Scalar k_tension,
                          Scalar k_bending);

#endif // SIMULABLE_GENERATOR_H_
