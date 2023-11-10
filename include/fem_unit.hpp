#ifndef FEM_UNIT_H_
#define FEM_UNIT_H_

#include "linear_algebra.hpp"
#include "physics_state.hpp"

struct FEM_UnitParameters {
    Scalar k, L0;
};

struct FEM_Unit {
    FEM_Unit(unsigned int p1,unsigned int p2, unsigned int p3, unsigned int p4, FEM_UnitParameters param);

    unsigned int p1, p2, p3, p4;
    FEM_UnitParameters parameters;

    Scalar get_energy(const PhysicsState& state);
    Vec3 get_force(const PhysicsState& state);
    Mat3 get_df_dx(const PhysicsState& state);
};

#endif // FEM_UNIT_H_
