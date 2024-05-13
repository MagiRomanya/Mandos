#ifndef SIMULABLE_GENERATOR_H_
#define SIMULABLE_GENERATOR_H_

#include "simulation.hpp"
#include "linear_algebra.hpp"

/**
 * Struct which stores indices from different data fields in the Simulation instance.
 *
 * For this to work, the simulable has to have contiguous degrees of freedom, particles and rigid bodies in the respective vectors.
 */
struct SimulableBounds {
    unsigned int dof_index, nDoF; // Starting dof index and total degrees of freedom of the simulable.
    unsigned int particle_index, n_particles; // Starting particle index and number of particles in the simulable.
    unsigned int rb_index, n_rb; // Starting rigid body index and number of rigid bodies in the simulable.
};

/**
 * Initializes a mass-spring simulable from a triangle mesh.
 *
 * The triangle mesh is interpreted as the rest configuration of the mass spring.
 * Each vertex of the mesh is interpreted as a particle and each vertex is interpreted as a tension spring.
 * This also creates extra springs to account for bending resistance, joining vertices separated by one edge.
 *
 * @param simulation The simulation instance where we want to create the simulable.
 * @param vertices, indices Description of the triangle indexed mesh.
 * @param node_mass Mass of each node.
 * @param k_tension, k_bending, damping Parameters related to elasticity and spring damping.
 */
SimulableBounds generate_mass_spring(Simulation& simulation,
                                     const std::vector<Scalar>& vertices,
                                     const std::vector<unsigned int>& indices,
                                     Scalar node_mass,
                                     Scalar k_tension,
                                     Scalar k_bending,
                                     Scalar damping);

/**
 * Initializes a single tetrahedron FEM soft body.
 *
 * This is a hard coded example and it is not possible to directly change the positions and velocities nor its rest configuration.
 *
 * @tparam MaterialType Type of FEM material.
 * @param simulation The simulation instance where we want to create the simulable.
 * @param node_mass Mass of each particle.
 * @param poisson_ratio, young_modulus Elasticity parameters of the soft body.
 */
template <typename MaterialType>
SimulableBounds generate_FEM3D_tetrahedron(Simulation& simulation, Scalar TotalMass, Scalar poisson_ratio, Scalar young_modulus);

/**
 * Initializes a FEM soft body from a tetrahedron mesh.
 *
 * The tetrahedron mesh is interpreted as the rest position of the soft body.
 * Each vertex is translated to a particle and each tetrahedron is translated into a FEM energy unit.
 *
 * @tparam MaterialType Type of FEM material.
 * @param simulation The simulation instance where we want to create the simulable.
 * @param node_mass Mass of each particle.
 * @param poisson_ratio, young_modulus Elasticity parameters of the soft body.
 * @param tet_indices, tet_vertices Description of the tetrahedron indexed mesh.
 */
template <typename MaterialType>
SimulableBounds generate_FEM3D_from_tetrahedron_mesh(Simulation& simulation, Scalar TotalMass, Scalar poisson_ratio, Scalar young_modulus,
                                                     const std::vector<unsigned int>& tet_indices, const std::vector<Scalar>& tet_vertices);


/**
 * Initializes a soft Rod from an origin and a given direction.
 *
 * @param simulation The simulation instance where we want to create the simulable.
 * @param segments The discretization of the rod. How many subdivisions should it have.
 * @param mass The total mass of the rod.
 * @param length The total length of the rod.
 * @param origin The starting point where the rod will be created.
 * @param direction Direction that the rod will grow from starting at the origin.
 * @param rod_parameters Physical parameters of the rod. (Note that L0 will be overwritten)
 */
SimulableBounds generate_rod(Simulation& simulation,
                             unsigned int segments, Scalar mass,
                             Scalar length, const Vec3 origin, const Vec3 direction,
                             const RodSegmentParameters& rod_parameters);

SimulableBounds generate_rod(Simulation& simulation,
                             std::vector<Scalar> vertices,
                             Scalar TotalMass,
                             const RodSegmentParameters& rod_parameters);


Vec3 compute_axis_angle_from_direction(const Vec3& direction);
#endif // SIMULABLE_GENERATOR_H_
