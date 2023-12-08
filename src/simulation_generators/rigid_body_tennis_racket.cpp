#include "inertia_energies.hpp"
#include "particle.hpp"
#include "raylib.h"
#include "rigid_body.hpp"
#include "simulable_generator.hpp"

SimulableBounds generate_RigidBody_tennis_racket_effect(Simulation& simulation) {
  const unsigned int nDoF = 6;
  // Mesh -> vertices & indices
  const Mesh RB_mesh = GenMeshPlane(8.0, 1.0, 20, 20);
  const unsigned int n_vertices = RB_mesh.vertexCount * 3;
  const unsigned int n_indices = RB_mesh.triangleCount * 3;
  std::vector<Scalar> vertices(RB_mesh.vertices, RB_mesh.vertices + n_vertices);
  std::vector<unsigned int> indices(RB_mesh.indices, RB_mesh.indices + n_indices);
  UnloadMesh(RB_mesh);

  const Scalar RB_MASS = 1.0;
  const Mat3 inertia_tensor = compute_initial_inertia_tensor_PARTICLES(RB_MASS, vertices);
  const unsigned int index = simulation.initial_state.x.size();
  simulation.initial_state.add_size(6);
  RigidBody rb(index, RB_MASS, inertia_tensor);
  add_rigid_body_to_simulation(simulation, rb);

  // Initial conditions
  simulation.initial_state.x.segment<nDoF>(index).setZero();
  simulation.initial_state.x_old.segment<nDoF>(index).setZero();
  set_angular_velocity(simulation.initial_state, simulation.TimeStep, rb.index+3, Vec3(0,0.01,0.1));

  return SimulableBounds{index, nDoF};
}
