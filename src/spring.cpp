#include "spring.hpp"
#include "linear_algebra.hpp"

Scalar ParticleSpring::compute_energy(Scalar TimeStep, const PhysicsState& state) const {
    const Vec3 x1 = p1.get_position(state);
    const Vec3 x2 = p2.get_position(state);
    const Scalar L = (x1 - x2).norm();

    const Scalar energy = parameters.get_energy(L);
    return  energy;
}

void ParticleSpring::compute_energy_and_derivatives(Scalar TimeStep, const PhysicsState& state, EnergyAndDerivatives& out) const {
    // Get the relevant sate
    // ---------------------------------------------------------------
    const Vec3 x1 = p1.get_position(state);
    const Vec3 x2 = p2.get_position(state);
    const Vec3 v1 = p1.get_velocity(state);
    const Vec3 v2 = p2.get_velocity(state);
    const Scalar L = (x1 - x2).norm();

    // Compute the energy derivatives
    // ---------------------------------------------------------------
    const Scalar energy = parameters.get_energy(L);
    const Vec3 gradient = parameters.get_energy_gradient(x1, x2, v1, v2, L);
    const Mat3 hessian = parameters.get_energy_hessian(TimeStep, x1, x2, L);

    // Add the energy derivatives to the global structure
    // ---------------------------------------------------------------
    out.energy += energy;
    for (unsigned int i = 0; i<3; i++) {
        // Newton's third law: equal and opposite reaction
        // Also force = - Grad E & df_dx = - Hess E
        out.gradient[p1.index + i] += gradient(i);
        out.gradient[p2.index + i] += -gradient(i);
        for (unsigned int j = 0; j<3; j++) {
            out.hessian_triplets.push_back(Triplet(p1.index+i, p1.index+j, hessian(i, j)));
            out.hessian_triplets.push_back(Triplet(p1.index+i, p2.index+j, -hessian(i, j)));
            out.hessian_triplets.push_back(Triplet(p2.index+i, p1.index+j, -hessian(i, j)));
            out.hessian_triplets.push_back(Triplet(p2.index+i, p2.index+j, hessian(i, j)));
        }
    }
}

Scalar SpringParameters::get_energy(Scalar L) const {
    return 0.5 * k * (L - L0) * (L - L0);
}

Vec3 SpringParameters::get_energy_gradient(const Vec3& x1, const Vec3& x2, const Vec3& v1, const Vec3& v2, Scalar L) const {
    /* Computes the spring force */
    const Vec3 u = (x1 - x2) / L;
    Vec3 f = -k * (L - L0) * u;
    // damping force
    f += - damping * u * u.transpose() * (v1 - v2);
    // The gradient is minus the force
    return -f;
}

Mat3 SpringParameters::get_energy_hessian(Scalar TimeStep, const Vec3& x1, const Vec3& x2, Scalar L) const {
    // u is the normalized vector between particles 1 and 2
    const Vec3 u = (x1 - x2) / L;
    // Initialize the derivative matrix
    Mat3 df_dx = (L - L0) * Mat3::Identity();
    const Mat3 uut = u * u.transpose();

    // Calculate the final derivative matrix
    df_dx = - k / L * (df_dx + L0 * uut);

    // Damping jacobian
    df_dx += - 1.0f / TimeStep * damping * uut;

    // The hessian is minus the force jacobian
    return -df_dx; // 3x3 matrix
}
