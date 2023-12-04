#ifndef SIMULATION_VISUALIZATION_H_
#define SIMULATION_VISUALIZATION_H_

#include <imgui.h>
#include "physics_render.hpp"
#include "physics_state.hpp"
#include "simulation.hpp"

Camera3D create_camera(unsigned int FPS = 120);

void simulation_render_simulables_and_energies(const Simulation& simulation, const PhysicsState& state);
void simulation_render_simulables(const Simulables& simulables, const PhysicsState& state);
void simulation_render_energies(const Energies& energies, const PhysicsState& state);

void simulation_visualization_loop(Simulation& simulation, PhysicsRenderers& phy_renderers);

#endif // SIMULATION_VISUALIZATION_H_
