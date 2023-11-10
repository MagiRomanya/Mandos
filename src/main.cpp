#include <iostream>
#include <memory>
#include <vector>

#include "linear_algebra.hpp"
#include "physics_state.hpp"
#include "spring.hpp"
#include "rigid_body.hpp"
#include "simulation.hpp"


int main(int argc, char *argv[]) {
    Vec a;
    PhysicsState state;
    state.resize(10);
    state.resize(10);
    std::cout << "Hello world!" << state.x.size() << std::endl;
    return 0;
}
