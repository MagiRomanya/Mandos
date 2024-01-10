#include "mandos.hpp"
#include "utility_functions.hpp"
#include "gravity.hpp"
#include "particle_rigid_body_copuling.hpp"
#include "spring.hpp"

RigidBodyHandle::RigidBodyHandle(Simulation& simulation, Scalar mass, const std::vector<Scalar> vertices)
    : rb(simulation.initial_state.get_nDoF(), mass, compute_initial_inertia_tensor_PARTICLES(mass, vertices)),
      rb_index(simulation.simulables.rigid_bodies.size()), simulation(simulation)
{
    simulation.initial_state.add_size(6);
    add_rigid_body_to_simulation(simulation, rb);
    simulation.initial_state.x.segment<6>(rb.index) = Eigen::Vector<Scalar,6>::Zero();
    simulation.initial_state.x_old.segment<6>(rb.index) = Eigen::Vector<Scalar,6>::Zero();
}

RigidBodyHandle::RigidBodyHandle(Simulation& simulation, Scalar mass, const Mat3& inertia_tensor)
    : rb(simulation.initial_state.get_nDoF(), mass, inertia_tensor),
      rb_index(simulation.simulables.rigid_bodies.size()), simulation(simulation)
{
    simulation.initial_state.add_size(6);
    add_rigid_body_to_simulation(simulation, rb);
    simulation.initial_state.x.segment<6>(rb.index) = Eigen::Vector<Scalar,6>::Zero();
    simulation.initial_state.x_old.segment<6>(rb.index) = Eigen::Vector<Scalar,6>::Zero();
}

RigidBodyHandle RigidBodyHandle::set_COM_initial_position(Vec3 pos) const {
    const Vec3 x = rb.get_COM_position(simulation.initial_state.x);
    const Vec3 x_old = rb.get_COM_position(simulation.initial_state.x_old);
    const Vec3 delta = x - x_old;
    simulation.initial_state.x.segment<3>(rb.index) = pos;
    simulation.initial_state.x_old.segment<3>(rb.index) = pos - delta;
    return *this;
}

RigidBodyHandle RigidBodyHandle::set_initial_orientation(Vec3 axis_angle) const {
    simulation.initial_state.x.segment<3>(rb.index + 3) = axis_angle;
    simulation.initial_state.x_old.segment<3>(rb.index + 3) = axis_angle;
    return *this;
}

RigidBodyHandle RigidBodyHandle::set_COM_initial_velocity(Vec3 vel) const {
    set_linear_velocity(simulation.initial_state, simulation.TimeStep, rb.index, vel);
    return *this;
}

RigidBodyHandle RigidBodyHandle::set_initial_angular_velocity(Vec3 omega) const {
    set_angular_velocity(simulation.initial_state, simulation.TimeStep, rb.index+3, omega);
    return *this;
}

RigidBodyHandle RigidBodyHandle::add_gravity(Scalar gravity) const {
    simulation.energies.gravities.emplace_back(rb.index+1, GravityParameters(gravity));
    return *this;
}

RigidBodyHandle RigidBodyHandle::freeze_translation() const {
    simulation.frozen_dof.push_back(rb.index+0);
    simulation.frozen_dof.push_back(rb.index+1);
    simulation.frozen_dof.push_back(rb.index+2);
    return *this;
}

RigidBodyHandle RigidBodyHandle::freeze_rotation() const {
    simulation.frozen_dof.push_back(rb.index+3+0);
    simulation.frozen_dof.push_back(rb.index+3+1);
    simulation.frozen_dof.push_back(rb.index+3+2);
    return *this;

}

Mat4 RigidBodyHandle::get_transformation_matrix(const PhysicsState& state) const {
    Eigen::Transform<Scalar, 3, Eigen::Affine> transformation;
    transformation.linear() = rb.compute_rotation_matrix(state.x);
    transformation.translation() = rb.get_COM_position(state.x);
    return transformation.matrix();
}

#ifdef ENABLE_LAGRANGE_MULTIPLIER_CONSTRAINTS
void join_rigid_bodies(Simulation& simulation, RigidBodyHandle rbA, Vec3 pA, RigidBodyHandle rbB, Vec3 pB) {
    simulation.constraints.rb_point_constraints.emplace_back(count_constraints(simulation.constraints), rbA.rb, rbB.rb, pA, pB);
}
#endif //ENABLE_LAGRANGE_MULTIPLIER_CONSTRAINTS

MassSpringHandle::MassSpringHandle(Simulation& simulation,
                                   const std::vector<Scalar>& vertices,
                                   const std::vector<unsigned int>& indices,
                                   Scalar TotalMass, Scalar k_tension, Scalar k_bending, Scalar damping)
    : TotalMass(TotalMass),
      bounds(generate_mass_spring(simulation, vertices, indices, 3*TotalMass / vertices.size(), k_tension, k_bending, damping)),
      simulation(simulation)
{}

Vec3 MassSpringHandle::compute_center_of_mass(const PhysicsState& state) const {
    Vec3 center_of_mass = Vec3::Zero();
    Scalar total_mass = 0;
    for (unsigned int i = bounds.particle_index; i < bounds.particle_index + bounds.n_particles; i++) {
        const Particle& p = simulation.simulables.particles[i];
        center_of_mass += p.get_position(state.x) * p.mass;
        total_mass += p.mass;
    }
    return center_of_mass / total_mass;
}

MassSpringHandle MassSpringHandle::set_initial_COM_position(const Vec3& position) const {
    const Vec3& COM = compute_center_of_mass(simulation.initial_state);
    const Vec3 displacement = position - COM;
    PhysicsState& state = simulation.initial_state;
    for (unsigned int i = bounds.particle_index; i < bounds.n_particles; i++) {
        const Particle& p = simulation.simulables.particles[i];
        state.x.segment<3>(p.index) = p.get_position(state.x) + displacement;
        state.x_old.segment<3>(p.index) = p.get_position(state.x_old) + displacement;
    }
    return *this;
}

MassSpringHandle MassSpringHandle::freeze_particles(const std::vector<unsigned int>& particle_indices) const {
    for (unsigned int i = 0; i < particle_indices.size(); i++) {
        const Particle& p = simulation.simulables.particles[i];
        simulation.frozen_dof.push_back(p.index+0);
        simulation.frozen_dof.push_back(p.index+1);
        simulation.frozen_dof.push_back(p.index+2);
    }
    return *this;
}

void MassSpringHandle::get_dof_vector(const PhysicsState& state, std::vector<float>& out_dofs) const {
    assert(out_dofs.size() == bounds.nDoF);
    for (unsigned int i = 0; i < bounds.nDoF; i++) {
        out_dofs[i] = state.x[i + bounds.dof_index];
    }
}

ParticleHandle::ParticleHandle(Simulation& simulation, Scalar mass)
    : particle(mass, simulation.initial_state.get_nDoF()),
      particle_index(simulation.simulables.particles.size()),
      simulation(simulation)
{
    simulation.initial_state.add_size(3);
    // Initialize initial conditions
    simulation.initial_state.x.segment<3>(particle.index) = Vec3::Zero();
    simulation.initial_state.x_old.segment<3>(particle.index) = Vec3::Zero();
    add_particle_to_simulation(simulation, particle);
}

ParticleHandle ParticleHandle::set_initial_position(Vec3 position) const {
    const Vec3 x = particle.get_position(simulation.initial_state.x);
    const Vec3 x_old = particle.get_position(simulation.initial_state.x_old);
    const Vec3 delta = x - x_old;
    simulation.initial_state.x.segment<3>(particle.index) = position;
    simulation.initial_state.x_old.segment<3>(particle.index) = position - delta;
    return *this;
}

ParticleHandle ParticleHandle::set_initial_velocity(Vec3 velocity) const {
    const Vec3 x = particle.get_position(simulation.initial_state.x);
    simulation.initial_state.x_old.segment<3>(particle.index) = x - velocity / simulation.TimeStep;
    return *this;
}

ParticleHandle ParticleHandle::add_gravity(Scalar gravity) const {
    simulation.energies.gravities.push_back(Gravity(particle.index+1, GravityParameters(gravity)));
    return *this;
}


ParticleHandle ParticleHandle::freeze() const {
    simulation.frozen_dof.push_back(particle.index+0);
    simulation.frozen_dof.push_back(particle.index+1);
    simulation.frozen_dof.push_back(particle.index+2);
    return *this;
}

void join_particles_with_spring(Simulation& simulation, const ParticleHandle& p1, const ParticleHandle& p2, Scalar k, Scalar damping) {
    const Vec3 x1 = p1.particle.get_position(simulation.initial_state.x);
    const Vec3 x2 = p2.particle.get_position(simulation.initial_state.x);
    const Scalar distance = (x1 -x2).norm();
    simulation.energies.particle_springs.emplace_back(p1.particle, p2.particle,
                                                      SpringParameters{.k = k, .L0 = distance, .damping = damping});
}

FEMHandle::FEMHandle(Simulation& simulation,
                     const std::vector<Scalar>& tetrahedron_vertices,
                     const std::vector<unsigned int>& tetrahedron_indices,
                     Scalar TotalMass, Scalar poisson_ratio, Scalar young_modulus)
    : TotalMass(TotalMass),
      bounds(generate_FEM3D_from_tetrahedron_mesh(simulation, 3 *TotalMass / tetrahedron_vertices.size(), poisson_ratio, young_modulus, tetrahedron_indices, tetrahedron_vertices)),
      simulation(simulation)
{}

Vec3 FEMHandle::compute_center_of_mass(const PhysicsState& state) const {
    Vec3 center_of_mass = Vec3::Zero();
    Scalar total_mass = 0;
    for (unsigned int i = bounds.particle_index; i < bounds.particle_index + bounds.n_particles; i++) {
        const Particle& p = simulation.simulables.particles[i];
        center_of_mass += p.get_position(state.x) * p.mass;
        total_mass += p.mass;
    }
    return center_of_mass / total_mass;
}

FEMHandle FEMHandle::freeze_particles(const std::vector<unsigned int>& particle_indices) const {
    for (unsigned int i = 0; i < particle_indices.size(); i++) {
        unsigned int index = particle_indices[i];
        simulation.frozen_dof.push_back(3*index+0 + bounds.dof_index);
        simulation.frozen_dof.push_back(3*index+1 + bounds.dof_index);
        simulation.frozen_dof.push_back(3*index+2 + bounds.dof_index);
    }
    return *this;
}


FEMHandle FEMHandle::add_gravity(Scalar gravity) const {
    for (unsigned int i = 0; i <bounds.nDoF/3; i++) {
        simulation.energies.gravities.push_back(Gravity(bounds.dof_index+3*i+1, GravityParameters(gravity)));
    }
    return *this;
}

void join_rigid_body_with_particle(Simulation& sim, RigidBodyHandle rb, ParticleHandle p) {
    ParticleRigidBodyCopuling copuling = ParticleRigidBodyCopuling(rb.rb, p.particle, p.particle.get_position(sim.initial_state.x));
    sim.copulings.add_copuling(copuling);
}
