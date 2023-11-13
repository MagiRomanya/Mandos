#ifndef EDIT_MODE_H_
#define EDIT_MODE_H_

#include <string>
#include <vector>

void edit_mode_render_meshes();

void edit_mode_sidebar();

struct MassSpringGUIGenerator {
    std::string mesh_name;
    float mass, k_tension, k_bending;
    std::vector<unsigned int> frozen_nodes;
};

#endif // EDIT_MODE_H_
