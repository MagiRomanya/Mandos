#include "spring.hpp"
#include "linear_algebra.hpp"

void ParticleSpring::compute_energy_and_derivatives(Scalar TimeStep, const PhysicsState& state, EnergyAndDerivatives& out) const {
    // Get the relevant sate
    // ---------------------------------------------------------------
    const Vec3 x1 = p1.get_position(state);
    const Vec3 x2 = p2.get_position(state);
    const Vec3 v1 = p1.get_velocity(state, TimeStep);
    const Vec3 v2 = p2.get_velocity(state, TimeStep);
    const Scalar L = (x1 - x2).norm();

    // Compute the energy derivatives
    // ---------------------------------------------------------------
    const Scalar energy = parameters.get_energy(L);
    const Vec3 force = parameters.get_force(x1, x2, v1, v2, L);
    const Mat3 df_dx = parameters.get_df_dx(TimeStep, x1, x2, v1, v2, L);

    // Add the energy derivatives to the global structure
    // ---------------------------------------------------------------
    out.energy += energy;
    for (unsigned int i = 0; i<3; i++) {
        out.jacobian[p1.index + i] += force(i);
        out.jacobian[p2.index + i] += -force(i);
        for (unsigned int j = 0; j<3; j++) {
            out.hessian_triplets.push_back(Triplet(p1.index+i, p1.index+j, df_dx(i, j)));
            out.hessian_triplets.push_back(Triplet(p1.index+i, p2.index+j, -df_dx(i, j)));
            out.hessian_triplets.push_back(Triplet(p2.index+i, p1.index+j, -df_dx(i, j)));
            out.hessian_triplets.push_back(Triplet(p2.index+i, p2.index+j, df_dx(i, j)));
        }
    }
}

Scalar SpringParameters::get_energy(Scalar L) const {
    return 0.5 * k * (L - L0) * (L - L0);
}

Vec3 SpringParameters::get_force(const Vec3& x1, const Vec3& x2, const Vec3& v1, const Vec3& v2, Scalar L) const {
    /* Computes the spring force */
    const Vec3 u = (x1 - x2) / L;
    Vec3 f = -k * (L - L0) * u;
    // damping force
    f += - damping * u * u.transpose() * (v1 - v2);
    return f;
}

Mat3 SpringParameters::get_df_dx(Scalar TimeStep, const Vec3& x1, const Vec3& x2, const Vec3& v1, const Vec3& v2, Scalar L) const {
    // u is the normalized vector between particles 1 and 2
    const Vec3 u = (x1 - x2) / L;
    // Initialize the derivative matrix
    Mat3 df_dx = (L - L0) * Mat3::Identity();
    // The u · u.transpose() matrix
    const Mat3 uut = u * u.transpose();
    // Calculate the final derivative matrix
    df_dx = - k / L * (df_dx + L0 * uut);

    // Damping jacobian
    df_dx += - 1/TimeStep * damping * uut;
    return df_dx; // 3x3 matrix
}
