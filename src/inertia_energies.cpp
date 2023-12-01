#include "inertia_energies.hpp"

void LinearInertia::compute_energy_and_derivatives(Scalar TimeStep, const PhysicsState& state, EnergyAndDerivatives& f) const {
    // Get the relevant sate
    // ---------------------------------------------------------------
    const Vec3 x = p.get_position(state);
    const Vec3 x_old = p.get_position_old(state);
    const Vec3 v_old = p.get_velocity_old(state, TimeStep);
    const Vec3 x_hat = x_old - TimeStep * v_old;
    const Scalar h2 = TimeStep*TimeStep;

    // Compute the energy derivatives
    // ---------------------------------------------------------------
    const Scalar energy = 1.0 / (2.0*h2) * (x - x_hat).transpose() * Mass * (x - x_hat);
    const Vec3 gradient = 1.0 / h2 * Mass * (x - x_hat);
    const Mat3 hessian = 1.0 / h2 * Mass;

    // Add the energy derivatives to the global structure
    // ---------------------------------------------------------------
    f.energy += energy;
    f.gradient.segment<3>(p.index) += gradient;
    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            f.hessian_triplets.emplace_back(p.index + i, p.index + j, hessian(i, j));
        }
    }
}

void RotationalInertia::compute_energy_and_derivatives(Scalar TimeStep, const PhysicsState& state, EnergyAndDerivatives& f) const {
    // Get the relevant sate
    // ---------------------------------------------------------------
    const Mat3 inertia_tensor = rb.inertia_tensor0;
    const Mat3 R = rb.compute_rotation_matrix(state.x);
    const Mat3 R_old = rb.compute_rotation_matrix(state.x_old);
    const Mat3 R_old2 = rb.compute_rotation_matrix(state.x_old2);
    const Mat3 R_guess = (R_old + R_old - R_old2);
    const Mat3 deltaR = R - R_guess;
    const Mat3 mat = R_guess * inertia_tensor * R.transpose();
    const Mat3 Asim = (mat - mat.transpose());
    const Mat3 Sim = (mat + mat.transpose()) / 2;
    const Scalar h2 = TimeStep * TimeStep;

    // Compute the energy derivatives
    // ---------------------------------------------------------------
    const Scalar KE = (deltaR * rb.inertia_tensor0 * deltaR.transpose()).trace();
    const Vec3 gradient = Vec3(Asim(1,2), - Asim(0,2), Asim(0,1)) / h2;
    const Mat3 hessian = 1.0 / h2 * (Sim.trace() * Mat3::Identity() - Sim);

    // Add the energy derivatives to the global structure
    // ---------------------------------------------------------------
    f.energy += KE;
    f.gradient.segment<3>(rb.index +3) += gradient;
    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            f.hessian_triplets.emplace_back(rb.index + 3 + i, rb.index + 3 + j, hessian(i, j));
        }
    }
}
