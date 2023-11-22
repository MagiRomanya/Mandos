#include "simulable_generator.hpp"
#include "simulation.hpp"

SimulableBounds generate_FEM3D_tetrahedron(Simulation& simulation, Scalar node_mass, Scalar poisson_ratio, Scalar young_modulus) {
  // Add 3 dimensions per particle to the degrees of freedom
  const unsigned int index = simulation.initial_state.x.size();
  simulation.initial_state.add_size(4*3);

  // Create particle simulables
  simulation.simulables.particles.emplace_back(node_mass, index+3*0);
  simulation.simulables.particles.emplace_back(node_mass, index+3*1);
  simulation.simulables.particles.emplace_back(node_mass, index+3*2);
  simulation.simulables.particles.emplace_back(node_mass, index+3*3);

  // Initial conditions
  simulation.initial_state.x.segment<3>(3*0+index) = Vec3(0,0,0);
  simulation.initial_state.x.segment<3>(3*1+index) = Vec3(1,0,0);
  simulation.initial_state.x.segment<3>(3*2+index) = Vec3(0,1,0);
  simulation.initial_state.x.segment<3>(3*3+index) = Vec3(0,0,1);
  simulation.initial_state.v.setZero();

  const Scalar mu = young_modulus / (2*(1 + poisson_ratio));
  const Scalar lambda = young_modulus * poisson_ratio / ((1+poisson_ratio) * (1-2*poisson_ratio));
  const Eigen::Matrix<Scalar,4,3> ds_dx =  compute_shape_function_derivative(Vec3(0,0,0), Vec3(1,0,0), Vec3(0,1,0), Vec3(0,0,1));
  FEM_ElementParameters FEM_paramters(mu, lambda, ds_dx);
  simulation.energies.fem_elements_3d.emplace_back(simulation.simulables.particles[0],
                                                   simulation.simulables.particles[1],
                                                   simulation.simulables.particles[2],
                                                   simulation.simulables.particles[3],
                                                   FEM_paramters);
  // Displace original position
  simulation.initial_state.x.segment<3>(3*2+index) = Vec3(0,1.5,0);

  return SimulableBounds{index, 12};
}
