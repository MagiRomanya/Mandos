#include "edge.hpp"
#include "gravity.hpp"
#include "mesh_boundary.hpp"
#include "particle.hpp"
#include "spring.hpp"
#include "simulable_generator.hpp"
#include <vector>

Scalar distance(const std::vector<Scalar>& vertices, unsigned int i1, unsigned int i2) {
    const Vec3 p1 = Vec3(vertices[3*i1], vertices[3*i1+1], vertices[3*i1+2]);
    const Vec3 p2 = Vec3(vertices[3*i2], vertices[3*i2+1], vertices[3*i2+2]);
    return (p1 - p2).norm();
}

SimulableBounds generate_mass_spring(Simulation& simulation,
                                     const std::vector<Scalar>& vertices,
                                     const std::vector<unsigned int>& indices,
                                     Scalar node_mass, Scalar k_tension,
                                     Scalar k_bending, Scalar damping)
{
    assert(vertices.size() % 3 == 0);
    const unsigned int index = simulation.initial_state.x.size();
    const unsigned int n_dof = vertices.size();

    // Resize degrees of freedom
    simulation.initial_state.add_size(n_dof);

    // Generate the particle simulables
    const unsigned int particle_index = simulation.simulables.particles.size();
    const std::vector<Particle>& particles = simulation.simulables.particles;
    for (unsigned int i = 0; i < n_dof; i+=3) {
        simulation.simulables.particles.emplace_back(node_mass, index + i);
    }
    // Create the inertia forces
    for (unsigned int i = particle_index; i < n_dof / 3; i++) {
        simulation.energies.linear_inertias.emplace_back(simulation.simulables.particles[i]);
    }

    // Initial conditions ( v = 0)
    for (size_t i=0; i < n_dof; i++) {
        simulation.initial_state.x[index + i] = vertices[i];
        simulation.initial_state.x_old[index + i] = vertices[i];
    }

    // Set up the springs
    std::vector<Edge> internalEdges, externalEdges;
    mesh_boundary(vertices, indices, internalEdges, externalEdges);

    const unsigned int n_flex = internalEdges.size() / 2.0 + externalEdges.size();
    const unsigned int n_bend = internalEdges.size() / 2.0;
    simulation.energies.particle_springs.reserve(simulation.energies.particle_springs.size() + n_flex + n_bend);

    for (size_t i = 0; i < externalEdges.size(); i++) {
        Edge &e = externalEdges[i];
        Scalar L0 = distance(vertices, e.a, e.b);
        SpringParameters param = {k_tension, L0, damping};
        simulation.energies.particle_springs.push_back(ParticleSpring(particles[particle_index + e.a], particles[particle_index + e.b], param));
    }

    for (size_t i = 0; i < internalEdges.size(); i += 2) {
        Edge &e1 = internalEdges[i];
        Edge &e2 = internalEdges[i + 1];
        Scalar L0 = distance(vertices, e1.a, e1.b);
        // Normal spring
        SpringParameters param = {k_tension, L0, damping};
        simulation.energies.particle_springs.push_back(ParticleSpring(particles[particle_index + e1.a], particles[particle_index + e1.b], param));

        // Bend spring
        L0 = distance(vertices, e1.opposite, e2.opposite);
        param = {k_bending, L0, damping};
        simulation.energies.particle_springs.push_back(ParticleSpring(particles[particle_index + e1.opposite], particles[particle_index + e2.opposite], param));
    }

    // Add gravity
    for (size_t i = 0; i < vertices.size(); i+=3) {
        GravityParameters param= {.intensity = static_cast<Scalar>(- node_mass)};
        // Add gravity to the y component
        simulation.energies.gravities.push_back(Gravity(index+i+1, param));
    }

    return SimulableBounds{.dof_index = index, .nDoF = n_dof};
};