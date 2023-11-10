#ifndef SIMULABLE_H_
#define SIMULABLE_H_

#include <memory>
#include <vector>
#include <cassert>

#include "rigid_body.hpp"
#include "fem_unit.hpp"
#include "spring.hpp"

struct Energies {
    std::vector<Spring> springs;
    std::vector<RigidBody> rigid_bodies;
    std::vector<FEM_Unit> fem_units;
};

struct Simulation {
    Energies energies;

    PhysicsState initial_state;

    std::vector<Triplet> initial_mass_matrix_triplets;
};

#endif // SIMULABLE_H_
