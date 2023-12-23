#include "inertia_energies.hpp"
#include "mesh.hpp"
#include "particle.hpp"
#include "rigid_body.hpp"
#include "simulable_generator.hpp"

SimulableBounds generate_RigidBody_tennis_racket_effect(Simulation& simulation) {
  const unsigned int nDoF = 6;
  const unsigned int rb_index = simulation.simulables.rigid_bodies.size();
  // Mesh -> vertices & indices
  SimulationMesh RB_mesh = SimulationMesh("resources/obj/cone.obj");
  const Scalar RB_MASS = 1.0;
  const Mat3 inertia_tensor = compute_initial_inertia_tensor_PARTICLES(RB_MASS, RB_mesh.vertices);
  const unsigned int index = simulation.initial_state.x.size();
  simulation.initial_state.add_size(6);
  const RigidBody rb(index, RB_MASS, inertia_tensor);
  add_rigid_body_to_simulation(simulation, rb);

  // Initial conditions
  simulation.initial_state.x.segment<nDoF>(index).setZero();
  simulation.initial_state.x_old.segment<nDoF>(index).setZero();
  set_angular_velocity(simulation.initial_state, simulation.TimeStep, rb.index+3, Vec3(0,0.01,0.1));

  return SimulableBounds{index, nDoF, 0, 0, rb_index, nDoF / 6};
}
