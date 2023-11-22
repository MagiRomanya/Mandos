#include "raylib.h"
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
  const Mat3 inertia_tensor = compute_initial_inertia_tensor(RB_MASS, indices, vertices, PARTICLES);
  const unsigned int index = simulation.initial_state.x.size();
  simulation.initial_state.add_size(6);
  RigidBody rb(index, RB_MASS, inertia_tensor);
  simulation.simulables.rigid_bodies.push_back(rb);

  // Initial conditions
  simulation.initial_state.x.segment<nDoF>(index).setZero();
  simulation.initial_state.v.segment<nDoF>(index).setZero();
  simulation.initial_state.v(index + 5) = 1; // add y direction angular velocity
  simulation.initial_state.v(index + 4) = 0.001; // add z direction angular velocity

  return SimulableBounds{index, nDoF};
}
