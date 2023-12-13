#include "inertia_energies.hpp"
#include "particle.hpp"
#include "simulation.hpp"

void LinearInertia::compute_energy_and_derivatives(Scalar TimeStep, const PhysicsState& state, const PhysicsState& state0, EnergyAndDerivatives& f) const {
    // Get the relevant sate
    // ---------------------------------------------------------------
    const Vec3 x = p.get_position(state.x);
    const Vec3 x0 = p.get_position(state0.x);
    const Vec3 v0 = p.get_velocity(state0, TimeStep);
    const Vec3 x_guess = x0 + TimeStep * v0;
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
    const Vec3 theta = state.x.segment<3>(rb.index+3);
    const Mat3 inertia_tensor = rb.J_inertia_tensor0;
    const Mat3 R = rb.compute_rotation_matrix(state.x);
    const Mat3 R0 = rb.compute_rotation_matrix(state0.x);
    const Mat3 R0old = rb.compute_rotation_matrix(state0.x_old);
    const Mat3 R_guess = (R0 + (R0 - R0old)); // x0 + h* v0

    const Mat3 deltaR = R - R_guess;
    const Mat3 rot_inertia = R * inertia_tensor * R_guess.transpose();
    const Mat3 S = (rot_inertia + rot_inertia.transpose()) / 2;
    const Mat3 A = (rot_inertia - rot_inertia.transpose()) / 2;
    const Scalar h2 = TimeStep * TimeStep;

    // Compute the energy derivatives
    // ---------------------------------------------------------------
    const Scalar KE = (deltaR * rb.J_inertia_tensor0 * deltaR.transpose()).trace();
    const Vec3 gradient = 2 * Vec3(-A(1,2), A(0,2), -A(0,1)) / h2;                 // v s.t. A = skew(v)
    const Mat3 hessian = 1.0 / h2 * (S.trace() * Mat3::Identity() - S);

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

void add_particle_to_simulation(Simulation& simulation, const Particle& p) {
    simulation.simulables.particles.push_back(p);
    simulation.energies.linear_inertias.emplace_back(p);
}

void add_rigid_body_to_simulation(Simulation& simulation, const RigidBody& rb) {
    simulation.simulables.rigid_bodies.emplace_back(rb);
    simulation.energies.linear_inertias.emplace_back(Particle(rb.mass, rb.index));
    simulation.energies.rotational_inertias.emplace_back(rb);
}
