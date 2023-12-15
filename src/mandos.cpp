#include "mandos.hpp"


RigidBodyHandle::RigidBodyHandle(Simulation& simulation, Scalar mass, const std::vector<Scalar> vertices)
    : rb(simulation.initial_state.get_nDoF(), mass, compute_initial_inertia_tensor_PARTICLES(mass, vertices)),
      rb_index(simulation.simulables.rigid_bodies.size())
{
    simulation.initial_state.add_size(6);
    add_rigid_body_to_simulation(simulation, rb);
    set_COM_initial_position(simulation, Vec3::Zero());
    set_initial_orientation(simulation, Vec3::Zero());
}

RigidBodyHandle::RigidBodyHandle(Simulation& simulation, Scalar mass, const Mat3& inertia_tensor)
    : rb(simulation.initial_state.get_nDoF(), mass, inertia_tensor),
      rb_index(simulation.simulables.rigid_bodies.size())
{
    simulation.initial_state.add_size(6);
    add_rigid_body_to_simulation(simulation, rb);
    set_COM_initial_position(simulation, Vec3::Zero());
    set_initial_orientation(simulation, Vec3::Zero());
}

RigidBodyHandle RigidBodyHandle::set_COM_initial_position(Simulation& simulation, Vec3 pos) const {
    simulation.initial_state.x.segment<3>(rb.index) = pos;
    simulation.initial_state.x_old.segment<3>(rb.index) = pos;
    return *this;
}

RigidBodyHandle RigidBodyHandle::set_initial_orientation(Simulation& simulation, Vec3 axis_angle) const {
    simulation.initial_state.x.segment<3>(rb.index + 3) = axis_angle;
    simulation.initial_state.x_old.segment<3>(rb.index + 3) = axis_angle;
    return *this;
}

RigidBodyHandle RigidBodyHandle::set_COM_initial_velocity(Simulation& simulation, Vec3 vel) const {
    set_linear_velocity(simulation.initial_state, simulation.TimeStep, rb.index, vel);
    return *this;
}

RigidBodyHandle RigidBodyHandle::set_initial_angluar_velocity(Simulation& simulation, Vec3 omega) const {
    set_angular_velocity(simulation.initial_state, simulation.TimeStep, rb.index, omega);
    return *this;
}

RigidBodyHandle RigidBodyHandle::add_gravity(Simulation& simulation, Scalar gravity) const {
    simulation.energies.gravities.emplace_back(rb.index+1, GravityParameters(gravity));
    return *this;
}

RigidBodyHandle RigidBodyHandle::freeze_translation(Simulation& simulation) const {
    simulation.frozen_dof.push_back(rb.index+0);
    simulation.frozen_dof.push_back(rb.index+1);
    simulation.frozen_dof.push_back(rb.index+2);
    return *this;
}

RigidBodyHandle RigidBodyHandle::freeze_rotation(Simulation& simulation) const {
    simulation.frozen_dof.push_back(rb.index+3+0);
    simulation.frozen_dof.push_back(rb.index+3+1);
    simulation.frozen_dof.push_back(rb.index+3+2);
    return *this;

}

Mat4 RigidBodyHandle::get_RB_transformation(const PhysicsState& state) const {
    Eigen::Transform<Scalar, 3, Eigen::Affine> transformation;
    transformation.linear() = rb.compute_rotation_matrix(state.x);
    transformation.translation() = rb.get_COM_position(state);
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
      bounds(generate_mass_spring(simulation, vertices, indices, 3*TotalMass / vertices.size(), k_tension, k_bending, damping))
{}

Vec3 MassSpringHandle::compute_center_of_mass(const Simulation& simulation, const PhysicsState& state) const {
    Vec3 center_of_mass = Vec3::Zero();
    Scalar total_mass = 0;
    for (unsigned int i = bounds.particle_index; i < bounds.particle_index + bounds.n_particles; i++) {
        const Particle& p = simulation.simulables.particles[i];
        center_of_mass += p.get_position(state.x) * p.mass;
        total_mass += p.mass;
    }
    return center_of_mass / total_mass;
}

MassSpringHandle MassSpringHandle::set_initial_COM_position(Simulation& simulation, const Vec3& position) const {
    const Vec3& COM = compute_center_of_mass(simulation, simulation.initial_state);
    const Vec3 displacement = position - COM;
    PhysicsState& state = simulation.initial_state;
    for (unsigned int i = bounds.particle_index; i < bounds.n_particles; i++) {
        const Particle& p = simulation.simulables.particles[i];
        state.x.segment<3>(p.index) = p.get_position(state.x) + displacement;
        state.x_old.segment<3>(p.index) = p.get_position(state.x_old) + displacement;
    }
    return *this;
}

MassSpringHandle MassSpringHandle::freeze_particles(Simulation& simulation, const std::vector<unsigned int>& particle_indices) const {
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
    : particle(simulation.initial_state.get_nDoF(), mass),
      particle_index(simulation.simulables.particles.size())
{
    simulation.initial_state.add_size(3);
    add_particle_to_simulation(simulation, particle);
}

ParticleHandle ParticleHandle::set_initial_position(Simulation& simulation, Vec3 position) const {
    const Vec3 x = particle.get_position(simulation.initial_state.x);
    const Vec3 x_old = particle.get_position(simulation.initial_state.x_old);
    const Vec3 delta = x - x_old;
    simulation.initial_state.x.segment<3>(particle.index) = position;
    simulation.initial_state.x.segment<3>(particle.index) = position - delta;
    return *this;
}

ParticleHandle ParticleHandle::set_initial_velocity(Simulation& simulation, Vec3 velocity) const {
    const Vec3 x = particle.get_position(simulation.initial_state.x);
    simulation.initial_state.x_old.segment<3>(particle.index) = x - velocity / simulation.TimeStep;
    return *this;
}


FEMHandle::FEMHandle(Simulation& simulation,
                     const std::vector<Scalar>& tetrahedron_vertices,
                     const std::vector<unsigned int>& tetrahedron_indices,
                     Scalar TotalMass, Scalar poisson_ratio, Scalar young_modulus)
    : TotalMass(TotalMass),
      bounds(generate_FEM3D_from_tetrahedron_mesh(simulation, 3 *TotalMass / tetrahedron_vertices.size(), poisson_ratio, young_modulus, tetrahedron_indices, tetrahedron_vertices))
{}

Vec3 FEMHandle::compute_center_of_mass(const Simulation& simulation, const PhysicsState& state) const {
    Vec3 center_of_mass = Vec3::Zero();
    Scalar total_mass = 0;
    for (unsigned int i = bounds.particle_index; i < bounds.particle_index + bounds.n_particles; i++) {
        const Particle& p = simulation.simulables.particles[i];
        center_of_mass += p.get_position(state.x) * p.mass;
        total_mass += p.mass;
    }
    return center_of_mass / total_mass;
}

FEMHandle FEMHandle::freeze_particles(Simulation& simulation, const std::vector<unsigned int>& particle_indices) const {
    for (unsigned int i = 0; i < particle_indices.size(); i++) {
        const Particle& p = simulation.simulables.particles[i];
        simulation.frozen_dof.push_back(p.index+0);
        simulation.frozen_dof.push_back(p.index+1);
        simulation.frozen_dof.push_back(p.index+2);
    }
    return *this;
}
