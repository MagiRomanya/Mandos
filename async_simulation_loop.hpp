#ifndef ASYNC_SIMULATION_LOOP_H_
#define ASYNC_SIMULATION_LOOP_H_

#include "physics_state.hpp"
#include "simulation.hpp"

void simulation_async_loop_request_iteration();

void simulation_async_loop(Simulation simulation);

void set_current_physics_state(PhysicsState state);

PhysicsState get_current_physics_state();

EnergyAndDerivatives get_current_energy_and_derivatives();

void simulation_async_end_thread();

#endif // ASYNC_SIMULATION_LOOP_H_
