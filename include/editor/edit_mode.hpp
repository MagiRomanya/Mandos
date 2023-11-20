#ifndef EDIT_MODE_H_
#define EDIT_MODE_H_

#include <string>
#include <vector>


enum GUI_STATE {
EDIT_MODE,
SIMULATION_MODE,
};

void edit_mode_loop();

struct MassSpringGUIGenerator {
    float mass, k_tension, k_bending;
    std::vector<unsigned int> frozen_nodes = {0,1,2};
};

struct PhysicsGUIGenerator {
    MassSpringGUIGenerator mass_spring;
};

#endif // EDIT_MODE_H_
