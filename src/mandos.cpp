#include <Eigen/Geometry>

#include "mandos.hpp"
#include "colliders.hpp"
#include "fem_element.hpp"
#include "inertia_energies.hpp"
#include "rigid_body.hpp"
#include "gravity.hpp"
#include "particle_rigid_body_copuling.hpp"
#include "rod_segment.hpp"
#include "simulable_generator.hpp"
#include "spring.hpp"
#include "utility_functions.hpp"

RigidBodyHandle::RigidBodyHandle(Simulation& simulation, Scalar mass, const std::vector<Scalar> vertices, bool global)
    : rb(simulation.initial_state.get_nDoF(), mass, compute_initial_inertia_tensor_PARTICLES(mass, vertices)),
      rb_index(static_cast<unsigned int>(simulation.simulables.rigid_bodies.size())),
      simulation(simulation)
{
    simulation.initial_state.add_size(6);
    add_rigid_body_to_simulation(simulation, rb, global);
    simulation.initial_state.x.segment<6>(rb.index) = Eigen::Vector<Scalar,6>::Zero();
    simulation.initial_state.v.segment<6>(rb.index) = Eigen::Vector<Scalar,6>::Zero();
}

RigidBodyHandle::RigidBodyHandle(Simulation& simulation, Scalar mass, const Mat3& inertia_tensor, bool global)
    : rb(simulation.initial_state.get_nDoF(), mass, inertia_tensor),
      rb_index(static_cast<unsigned int>(simulation.simulables.rigid_bodies.size())),
      simulation(simulation)
{
    simulation.initial_state.add_size(6);
    add_rigid_body_to_simulation(simulation, rb, global);
    simulation.initial_state.x.segment<6>(rb.index) = Eigen::Vector<Scalar,6>::Zero();
    simulation.initial_state.v.segment<6>(rb.index) = Eigen::Vector<Scalar,6>::Zero();
}

RigidBodyHandle RigidBodyHandle::set_COM_initial_position(Vec3 pos) const {
    const Vec3 x = rb.get_COM_position(simulation.initial_state.x);
    simulation.initial_state.x.segment<3>(rb.index) = pos;
    return *this;
}

RigidBodyHandle RigidBodyHandle::set_initial_orientation(Vec3 axis_angle) const {
    simulation.initial_state.x.segment<3>(rb.index + 3) = clamp_axis_angle(axis_angle);
    return *this;
}

RigidBodyHandle RigidBodyHandle::set_COM_initial_velocity(Vec3 vel) const {
    simulation.initial_state.v.segment<3>(rb.index) = vel;
    return *this;
}

RigidBodyHandle RigidBodyHandle::set_initial_angular_velocity(Vec3 omega) const {
    simulation.initial_state.v.segment<3>(rb.index+3) = omega;
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

Scalar RigidBodyHandle::compute_energy(const PhysicsState& state, const PhysicsState& state0) {
    RotationalInertia i = RotationalInertia(rb);
    return i.compute_energy(simulation.TimeStep, state, state0);
}

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
        center_of_mass += p.get_position(state) * p.mass;
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
        state.x.segment<3>(p.index) = p.get_position(state) + displacement;
    }
    return *this;
}

MassSpringHandle MassSpringHandle::freeze_particles(const std::vector<unsigned int>& particle_indices) const {
    for (unsigned int i = 0; i < particle_indices.size(); i++) {
        unsigned int index = particle_indices[i];
        simulation.frozen_dof.push_back(3*index+0 + bounds.dof_index);
        simulation.frozen_dof.push_back(3*index+1 + bounds.dof_index);
        simulation.frozen_dof.push_back(3*index+2 + bounds.dof_index);
    }
    return *this;
}

MassSpringHandle MassSpringHandle::add_gravity(Scalar gravity) const {
    for (unsigned int i = 0; i <bounds.nDoF/3; i++) {
        simulation.energies.gravities.push_back(Gravity(bounds.dof_index+3*i+1, GravityParameters(gravity)));
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
      particle_index(static_cast<unsigned int>(simulation.simulables.particles.size())),
      simulation(simulation)
{
    simulation.initial_state.add_size(3);
    // Initialize initial conditions
    simulation.initial_state.x.segment<3>(particle.index) = Vec3::Zero();
    simulation.initial_state.v.segment<3>(particle.index) = Vec3::Zero();
    add_particle_to_simulation(simulation, particle);
}

ParticleHandle ParticleHandle::set_initial_position(Vec3 position) const {
    const Vec3 x = particle.get_position(simulation.initial_state);
    simulation.initial_state.x.segment<3>(particle.index) = position;
    return *this;
}

ParticleHandle ParticleHandle::set_initial_velocity(Vec3 velocity) const {
    simulation.initial_state.v.segment<3>(particle.index) = velocity;
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
    const Vec3 x1 = p1.particle.get_position(simulation.initial_state);
    const Vec3 x2 = p2.particle.get_position(simulation.initial_state);
    const Scalar distance = (x1 -x2).norm();
    simulation.energies.particle_springs.emplace_back(p1.particle, p2.particle,
                                                      SpringParameters{.k = k, .L0 = distance, .damping = damping});
}

FEMHandle::FEMHandle(Simulation& simulation,
                     const std::vector<Scalar>& tetrahedron_vertices,
                     const std::vector<unsigned int>& tetrahedron_indices,
                     Scalar TotalMass, Scalar poisson_ratio, Scalar young_modulus)
    : TotalMass(TotalMass),
      bounds(generate_FEM3D_from_tetrahedron_mesh<FEM_NeoHookeanMaterial>(simulation, 3 *TotalMass / tetrahedron_vertices.size(), poisson_ratio, young_modulus, tetrahedron_indices, tetrahedron_vertices)),
      simulation(simulation)
{}

Vec3 FEMHandle::compute_center_of_mass(const PhysicsState& state) const {
    Vec3 center_of_mass = Vec3::Zero();
    Scalar total_mass = 0;
    for (unsigned int i = bounds.particle_index; i < bounds.particle_index + bounds.n_particles; i++) {
        const Particle& p = simulation.simulables.particles[i];
        center_of_mass += p.get_position(state) * p.mass;
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

ParticleHandle get_particle_handle(Simulation& sim, unsigned int particle_index) {
    Particle p = sim.simulables.particles[particle_index];
    return ParticleHandle(sim, p, particle_index);
}

void join_rigid_body_with_particle(Simulation& sim, RigidBodyHandle rb, ParticleHandle p) {
    ParticleRigidBodyCopuling copuling = ParticleRigidBodyCopuling(rb.rb, p.particle, p.particle.get_position(sim.initial_state));
    sim.copulings.add_copuling(copuling);
}

RodHandle::RodHandle(Simulation& simulation, unsigned int segments, const Vec3& rod_vector, Scalar TotalMass, const RodSegmentParameters& parameters)
    : simulation(simulation), TotalMass(TotalMass),
      bounds(generate_rod(simulation, segments, TotalMass, rod_vector.norm(), Vec3::Zero(), rod_vector.normalized(), parameters))
{}

RodHandle::RodHandle(Simulation& simulation, const std::vector<Scalar>& vertices, Scalar TotalMass, const RodSegmentParameters& parameters)
    : simulation(simulation), TotalMass(TotalMass),
      bounds(generate_rod(simulation, vertices, TotalMass, parameters))
{}

Vec3 RodHandle::compute_center_of_mass(const PhysicsState& state) const {
    Vec3 com = Vec3::Zero();
    for (unsigned int i = 0; i < bounds.n_rb; i++) {
        com += simulation.simulables.rigid_bodies[bounds.rb_index + i].get_COM_position(state.x);
    }
    com /= bounds.n_rb;
    return com;
}

RodHandle RodHandle::add_gravity(Scalar gravity) const {
    for (unsigned int i = 0; i < bounds.n_rb; i++) {
        const RigidBody& rb = simulation.simulables.rigid_bodies[bounds.rb_index + i];
        simulation.energies.gravities.emplace_back(rb.index + 1, GravityParameters(gravity));
    }

    return *this;
}

RodHandle RodHandle::set_initial_rod_position(const Vec3& origin) const {
    for (unsigned int i = 0; i < bounds.n_rb; i++) {
        const RigidBody& rb = simulation.simulables.rigid_bodies[bounds.rb_index + i];
        simulation.initial_state.x.segment<3>(rb.index) += origin;
    }
    return *this;
}

RodHandle RodHandle::set_initial_rigid_body_position(const Scalar s, const Vec3& x) const {
    assert(s >= 0 && s <= 1.0);
    const unsigned int index = static_cast<unsigned int>(s * (bounds.n_rb));
    const RigidBody& rb = simulation.simulables.rigid_bodies[bounds.rb_index + index];
    simulation.initial_state.x.segment<3>(rb.index) = x;
    return *this;
}

RodHandle RodHandle::set_initial_rigid_body_velocity(const Scalar s, const Vec3& v) const {
    assert(s >= 0 && s <= 1.0);
    const unsigned int index = static_cast<unsigned int>(s * (bounds.n_rb - 1));
    const RigidBody& rb = simulation.simulables.rigid_bodies[bounds.rb_index + index];
    simulation.initial_state.v.segment<3>(rb.index) = v;
    return *this;
}

// RodHandle RodHandle::set_initial_rod_direction(const Vec3& direction) const {
//     // Changing direction means to change the postions of the rigid bodies
//     // and also their orientation

//     // Compute the orientation:
//     Vec3 normalized_direction = direction.normalized();
//     Vec3 axis_angle = Vec3::Zero();
//     if (not normalized_direction.isApprox(Vec3(0.0, 0.0, 1.0), 1e-6)) {
//         const Vec3 tangent = cross(normalized_direction, Vec3(0.0, 1.0, 0.0)).normalized();
//         const Vec3 bitangent = cross(normalized_direction, tangent).normalized();
//         Mat3 rotation;
//         rotation << tangent, bitangent, normalized_direction;
//         Eigen::AngleAxis<Scalar> angle_axis = Eigen::AngleAxis<Scalar>(rotation);
//         axis_angle = angle_axis.axis() * angle_axis.angle();
//     }

//     const RigidBody& rbOrigin = simulation.simulables.rigid_bodies[bounds.rb_index];
//     const Vec3 origin = rbOrigin.get_COM_position(simulation.initial_state.x);

//     for (unsigned int i = 0; i < bounds.n_rb; i++) {
//         const RigidBody& rb = simulation.simulables.rigid_bodies[bounds.rb_index + i];
//         simulation.initial_state.x.segment<3>(rb.index) = origin + L0 * normalized_direction * i;
//         simulation.initial_state.x.segment<3>(rb.index+3) = axis_angle;
//     }

//     return *this;
// }

RodHandle RodHandle::freeze_rigid_body(Scalar s) const {
    assert(s >= 0 && s <= 1.0);
    const unsigned int index = static_cast<unsigned int>(s * (bounds.n_rb));
    const RigidBody& rb = simulation.simulables.rigid_bodies[bounds.rb_index + index];
    for (unsigned int i = 0; i < 6; i++)
        simulation.frozen_dof.push_back(rb.index + index + i);

    return *this;
}


PlaneColliderHandle::PlaneColliderHandle(Simulation& simulation)
    :simulation(simulation)
{
    PlaneCollider collider;
    collider.center = Vec3::Zero();
    collider.normal = Vec3(0.0, 1.0, 0.0);
    index = simulation.colliders.plane_colliders.size();
    simulation.colliders.plane_colliders.push_back(collider);
}

PlaneColliderHandle PlaneColliderHandle::set_origin_position(const Vec3& origin) const {
    simulation.colliders.plane_colliders[index].center = origin;
    return *this;
}

PlaneColliderHandle PlaneColliderHandle::set_direction(const Vec3& direction) const {
    simulation.colliders.plane_colliders[index].normal = direction.normalized();
    return *this;
}

SphereColliderHandle::SphereColliderHandle(Simulation& simulation)
    :simulation(simulation)
{
    SphereCollider collider;
    collider.center = Vec3::Zero();
    collider.radius = 10.0;
    index = simulation.colliders.sphere_colliders.size();
    simulation.colliders.sphere_colliders.push_back(collider);
}

SphereColliderHandle SphereColliderHandle::set_origin_position(const Vec3& origin) const {
    simulation.colliders.sphere_colliders[index].center = origin;
    return *this;
}

SphereColliderHandle SphereColliderHandle::set_radius(const Scalar& radius) const {
    simulation.colliders.sphere_colliders[index].radius = radius;
    return *this;
}


SDFColliderHandle::SDFColliderHandle(Simulation& simulation, const SimulationMesh& sim_mesh)
    :simulation(simulation)
{
    index = simulation.colliders.sdf_colliders.size();
    simulation.colliders.sdf_colliders.emplace_back(sim_mesh);
}
