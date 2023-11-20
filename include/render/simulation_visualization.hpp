#ifndef SIMULATION_VISUALIZATION_H_
#define SIMULATION_VISUALIZATION_H_

#include <imgui.h>
#include "physics_render.hpp"
#include "simulation.hpp"

Camera3D create_camera(unsigned int FPS = 60);

void simulation_visualization_loop(Simulation& simulation, PhysicsRenderers& phy_renderers);

#endif // SIMULATION_VISUALIZATION_H_
