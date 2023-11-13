#include "spring.hpp"
#include "linear_algebra.hpp"
#include <iostream>

void Spring::compute_energy_and_derivatives(const PhysicsState& state, EnergyAndDerivatives& out) const {
    Vec3 x1 = Vec3(state.x[p1], state.x[p1+1], state.x[p1+2]);
    Vec3 x2 = Vec3(state.x[p2], state.x[p2+1], state.x[p2+2]);
    Scalar L = (x1 - x2).norm();

    out.energy += get_energy(x1, x2, L);
    Vec3 force = get_force(x1, x2, L);
    Mat3 df_dx = get_df_dx(x1, x2, L);

    for (unsigned int i = 0; i<3; i++) {
        out.force[p1 + i] += force(i);
        out.force[p2 + i] += -force(i);
        for (unsigned int j = 0; j<3; j++) {
            out.df_dx_triplets.push_back(Triplet(p1+i, p1+j, df_dx(i, j)));
            out.df_dx_triplets.push_back(Triplet(p1+i, p2+j, -df_dx(i, j)));
            out.df_dx_triplets.push_back(Triplet(p2+i, p1+j, -df_dx(i, j)));
            out.df_dx_triplets.push_back(Triplet(p2+i, p2+j, df_dx(i, j)));
        }
    }
}

Scalar Spring::get_energy(const Vec3& x1, const Vec3& x2, Scalar L) const {
    return 0.5 * parameters.k * (L - parameters.L0) * (L - parameters.L0);
}

Vec3 Spring::get_force(const Vec3& x1, const Vec3& x2, Scalar L) const {
    /* Computes the spring force */
    Vec3 f = -parameters.k * (L - parameters.L0) * (x1 - x2) / L;
    return f;
}

Mat3 Spring::get_df_dx(const Vec3& x1, const Vec3& x2, Scalar L) const {
    // u is the normalized vector between particles 1 and 2
    Vec3 u = (x1 - x2) / L;
    // Initialize the derivative matrix
    Mat3 df_dx = (L - parameters.L0) * Mat3::Identity();
    // The u · u.transpose() matrix
    Mat3 uut = u * u.transpose();

    // Calculate the final derivative matrix
    df_dx = - parameters.k / L * (df_dx + parameters.L0 * uut);

    return df_dx; // 3x3 matrix
}
