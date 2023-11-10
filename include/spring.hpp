#ifndef SPRING_H_
#define SPRING_H_

#include "linear_algebra.hpp"
#include "physics_state.hpp"

struct SpringParameters {
    Scalar k, L0;
};

struct Spring {
    Spring(unsigned int p1, unsigned int p2, SpringParameters param)
        : p1(p1), p2(p2), parameters(param) {};

    unsigned int p1, p2;
    SpringParameters parameters;

    void compute_energy_and_derivatives(const PhysicsState& state, EnergyAndDerivatives& out) const;

    private:
        Scalar get_energy(const Vec3& x1, const Vec3& x2, Scalar L) const;
        Vec3 get_force(const Vec3& x1, const Vec3& x2, Scalar L) const;
        Mat3 get_df_dx(const Vec3& x1, const Vec3& x2, Scalar L) const;
};

#endif // SPRING_H_
