#include "inertia_energies.hpp"
#include "particle.hpp"
#include "simulable_generator.hpp"
#include "simulation.hpp"

SimulableBounds generate_FEM3D_tetrahedron(Simulation& simulation, Scalar node_mass, Scalar poisson_ratio, Scalar young_modulus) {
  // Add 3 dimensions per particle to the degrees of freedom
  const unsigned int index = simulation.initial_state.x.size();
  const unsigned int nDoF = 12;
  simulation.initial_state.add_size(nDoF);

  // Create particle simulables
  const unsigned int particle_index = simulation.simulables.particles.size();
  add_particle_to_simulation(simulation, Particle(node_mass, index + 3*0));
  add_particle_to_simulation(simulation, Particle(node_mass, index + 3*1));
  add_particle_to_simulation(simulation, Particle(node_mass, index + 3*2));
  add_particle_to_simulation(simulation, Particle(node_mass, index + 3*3));

  // Initial conditions
  simulation.initial_state.x.segment<3>(3*0+index) = Vec3(0,0,0);
  simulation.initial_state.x.segment<3>(3*1+index) = Vec3(1,0,0);
  simulation.initial_state.x.segment<3>(3*2+index) = Vec3(0,1,0);
  simulation.initial_state.x.segment<3>(3*3+index) = Vec3(0,0,1);

  simulation.initial_state.x_old.segment<nDoF>(index) = simulation.initial_state.x.segment<nDoF>(index);

  const Scalar mu = young_modulus / (2*(1 + poisson_ratio));
  const Scalar lambda = young_modulus * poisson_ratio / ((1+poisson_ratio) * (1-2*poisson_ratio));
  const Eigen::Matrix<Scalar,4,3> ds_dx =  compute_shape_function_derivative(Vec3(0,0,0), Vec3(1,0,0), Vec3(0,1,0), Vec3(0,0,1));
  FEM_ElementParameters FEM_paramters(mu, lambda, ds_dx);
  simulation.energies.fem_elements_3d.emplace_back(simulation.simulables.particles[particle_index + 0],
                                                   simulation.simulables.particles[particle_index + 1],
                                                   simulation.simulables.particles[particle_index + 2],
                                                   simulation.simulables.particles[particle_index + 3],
                                                   FEM_paramters);
  // Displace original position
  simulation.initial_state.x.segment<3>(3*2+index) = Vec3(0,1.5,0);
  simulation.initial_state.x_old.segment<3>(3*2+index) = Vec3(0,1.5,0);

  return SimulableBounds{index, nDoF};
}

SimulableBounds generate_FEM3D_from_tetrahedron_mesh(Simulation& simulation, Scalar node_mass, Scalar poisson_ratio, Scalar young_modulus,
                                                     const std::vector<unsigned int>& tet_indices, const std::vector<float>& tet_vertices) {
  const unsigned int index = simulation.initial_state.x.size();
  const unsigned int particle_index = simulation.simulables.particles.size();
  const unsigned int nDoF = tet_vertices.size();
  const unsigned int n_particles = nDoF / 3;
  const unsigned int n_tets = tet_indices.size() / 4;

  // Allocate space for the simualtion state
  simulation.initial_state.add_size(nDoF);

  // Create particle simulables
  for (unsigned int i = 0; i < n_particles; i++) {
    add_particle_to_simulation(simulation, Particle(node_mass, index + 3*i));
  }

  // Initial conditions
  for (unsigned int i = 0; i < nDoF; i++) {
    simulation.initial_state.x[index+i] = tet_vertices[i];
    simulation.initial_state.x_old[index+i] = tet_vertices[i];
  }

  const Scalar mu = young_modulus / (2*(1 + poisson_ratio));
  const Scalar lambda = young_modulus * poisson_ratio / ((1+poisson_ratio) * (1-2*poisson_ratio));

  // Generate FEM element for each tetrahedron
  for (unsigned int i=0; i < n_tets; i++) {
    const Particle p1 = simulation.simulables.particles[particle_index + tet_indices[4*i+0]];
    const Particle p2 = simulation.simulables.particles[particle_index + tet_indices[4*i+1]];
    const Particle p3 = simulation.simulables.particles[particle_index + tet_indices[4*i+2]];
    const Particle p4 = simulation.simulables.particles[particle_index + tet_indices[4*i+3]];
    const Eigen::Matrix<Scalar,4,3> ds_dx =  compute_shape_function_derivative(p1.get_position(simulation.initial_state.x),
                                                                               p2.get_position(simulation.initial_state.x),
                                                                               p3.get_position(simulation.initial_state.x),
                                                                               p4.get_position(simulation.initial_state.x));
    FEM_ElementParameters FEM_paramters(mu, lambda, ds_dx);
    simulation.energies.fem_elements_3d.emplace_back(p1, p2, p3, p4, FEM_paramters);
  }

  // Displace original position
  simulation.initial_state.x[index] += 2;
  simulation.initial_state.x_old[index] += 2;

  return SimulableBounds{index, nDoF};
}
