#ifndef SIMULABLE_GENERATOR_H_
#define SIMULABLE_GENERATOR_H_

#include "simulation.hpp"
#include "linear_algebra.hpp"

struct SimulableBounds {
    unsigned int dof_index, nDoF;
    unsigned int particle_index, n_particles;
    unsigned int rb_index, n_rb;
};

SimulableBounds generate_mass_spring(Simulation& simulation,
                                     const std::vector<Scalar>& vertices,
                                     const std::vector<unsigned int>& indices,
                                     Scalar node_mass,
                                     Scalar k_tension,
                                     Scalar k_bending,
                                     Scalar damping);

template <typename MaterialType>
SimulableBounds generate_FEM3D_tetrahedron(Simulation& simulation, Scalar node_mass, Scalar poisson_ratio, Scalar young_modulus);

template <typename MaterialType>
SimulableBounds generate_FEM3D_from_tetrahedron_mesh(Simulation& simulation, Scalar node_mass, Scalar poisson_ratio, Scalar young_modulus,
                                                     const std::vector<unsigned int>& tet_indices, const std::vector<Scalar>& tet_vertices);

SimulableBounds generate_RigidBody_tennis_racket_effect(Simulation& simulation);

#endif // SIMULABLE_GENERATOR_H_
