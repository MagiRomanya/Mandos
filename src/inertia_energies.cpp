#include "inertia_energies.hpp"

void LinearInertia::compute_energy_and_derivatives(Scalar TimeStep, const PhysicsState& state, const PhysicsState& state0, EnergyAndDerivatives& f) const {
    // Get the relevant sate
    // ---------------------------------------------------------------
    const Vec3 x = p.get_position(state.x);
    const Vec3 x0 = p.get_position(state0.x);
    const Vec3 v0 = p.get_velocity(state0, TimeStep);
    const Vec3 x_guess = x0 - TimeStep * v0;
    const Scalar h2 = TimeStep*TimeStep;

    // Compute the energy derivatives
    // ---------------------------------------------------------------
    const Scalar energy = 1.0 / (2.0*h2) * (x - x_guess).transpose() * Mass * (x - x_guess);
    const Vec3 gradient = 1.0 / h2 * Mass * (x - x_guess);
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

void RotationalInertia::compute_energy_and_derivatives(Scalar TimeStep, const PhysicsState& state, const PhysicsState& state0, EnergyAndDerivatives& f) const {
    // Get the relevant sate
    // ---------------------------------------------------------------
    const Mat3 inertia_tensor = rb.inertia_tensor0;
    const Mat3 R = rb.compute_rotation_matrix(state.x);
    const Mat3 R0 = rb.compute_rotation_matrix(state0.x);
    const Mat3 R0old = rb.compute_rotation_matrix(state0.x_old);
    const Mat3 R_guess = (R0 + (R0 - R0old)); // x0 + h* v0

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
