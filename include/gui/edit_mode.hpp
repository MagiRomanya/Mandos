#ifndef EDIT_MODE_H_
#define EDIT_MODE_H_

#include "raylib.h"
#include "simulable_generator.hpp"
#include "simulation.hpp"
#include <string>
#include <vector>

void edit_mode_loop();

void edit_mode_render_meshes();

void edit_mode_sidebar();


struct MassSpringGUIGenerator {
    float mass, k_tension, k_bending;
    std::vector<unsigned int> frozen_nodes = {0,1,2};
};

struct PhysicsGUIGenerator {
    MassSpringGUIGenerator mass_spring;
};

#endif // EDIT_MODE_H_
