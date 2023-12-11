#ifndef MANDOS_H_
#define MANDOS_H_

#include <vector>

#include "gravity.hpp"
#include "linear_algebra.hpp"
#include "inertia_energies.hpp"
#include "particle.hpp"
#include "physics_state.hpp"
#include "rigid_body.hpp"
#include "simulable_generator.hpp"
#include "simulation.hpp"
#include "spring.hpp"

class RigidBodyHandle {
    public:
        RigidBodyHandle(Simulation& simulation, Scalar mass, const std::vector<Scalar> vertices)
            : rb_index(simulation.simulables.rigid_bodies.size()),
              rb(simulation.initial_state.get_nDoF(), mass, compute_initial_inertia_tensor_PARTICLES(mass, vertices))
        {
            simulation.initial_state.add_size(6);
            add_rigid_body_to_simulation(simulation, rb);
            set_COM_initial_position(simulation, Vec3::Zero());
            set_initial_orientation(simulation, Vec3::Zero());
        }

        RigidBodyHandle(Simulation& simulation, Scalar mass, const Mat3& inertia_tensor)
            : rb_index(simulation.simulables.rigid_bodies.size()),
              rb(simulation.initial_state.get_nDoF(), mass, inertia_tensor)
        {
            simulation.initial_state.add_size(6);
            add_rigid_body_to_simulation(simulation, rb);
            set_COM_initial_position(simulation, Vec3::Zero());
            set_initial_orientation(simulation, Vec3::Zero());
        }

        inline RigidBodyHandle set_COM_initial_position(Simulation& simulation, Vec3 pos) const {
            simulation.initial_state.x.segment<3>(rb.index) = pos;
            simulation.initial_state.x_old.segment<3>(rb.index) = pos;
            return *this;
        }

        inline RigidBodyHandle set_initial_orientation(Simulation& simulation, Vec3 axis_angle) const {
            simulation.initial_state.x.segment<3>(rb.index + 3) = axis_angle;
            simulation.initial_state.x_old.segment<3>(rb.index + 3) = axis_angle;
            return *this;
        }

        inline RigidBodyHandle set_COM_initial_velocity(Simulation& simulation, Vec3 vel) const {
            set_linear_velocity(simulation.initial_state, simulation.TimeStep, rb.index, vel);
            return *this;
        }
        inline RigidBodyHandle set_initial_angluar_velocity(Simulation& simulation, Vec3 omega) const {
            set_angular_velocity(simulation.initial_state, simulation.TimeStep, rb.index, omega);
            return *this;
        }
        inline RigidBodyHandle add_gravity(Simulation& simulation, Scalar gravity) const {
            simulation.energies.gravities.emplace_back(rb.index+1, GravityParameters(gravity));
            return *this;
        }

        inline Mat4 get_RB_transformation(const PhysicsState& state) const {
            Eigen::Transform<Scalar, 3, Eigen::Affine> transformation;
            transformation.linear() = rb.compute_rotation_matrix(state.x);
            transformation.translation() = rb.get_COM_position(state);
            return transformation.matrix();
        }

        const RigidBody rb;
        const unsigned int rb_index;
};


class MassSpringHandle {
    public:
        MassSpringHandle(Simulation& simulation,
                         const std::vector<Scalar>& vertices,
                         const std::vector<unsigned int>& indices,
                         Scalar TotalMass, Scalar k_tension, Scalar k_bending, Scalar damping)
            : bounds(generate_mass_spring(simulation, vertices, indices, 3*TotalMass / vertices.size(), k_tension, k_bending, damping)),
              TotalMass(TotalMass)
        {}

        inline Vec3 compute_center_of_mass(const Simulation& simulation, const PhysicsState& state) const {
            Vec3 center_of_mass = Vec3::Zero();
            Scalar total_mass = 0;
            for (unsigned int i = bounds.particle_index; i < bounds.particle_index + bounds.n_particles; i++) {
                const Particle& p = simulation.simulables.particles[i];
                center_of_mass += p.get_position(state.x) * p.mass;
                total_mass += p.mass;
            }
            return center_of_mass / total_mass;
        }

        inline unsigned int get_n_particles() const { return bounds.n_particles; }

        inline void freeze_particles(Simulation& simulation, const std::vector<unsigned int>& particle_indices) const {
            for (unsigned int i = 0; i < particle_indices.size(); i++) {
                const Particle& p = simulation.simulables.particles[i];
                simulation.frozen_dof.push_back(p.index+0);
                simulation.frozen_dof.push_back(p.index+1);
                simulation.frozen_dof.push_back(p.index+2);
            }
        }

        const Scalar TotalMass;
        const SimulableBounds bounds;
};

class FEMHandle {
    public:
        FEMHandle(Simulation& simulation,
                  const std::vector<Scalar>& tetrahedron_vertices,
                  const std::vector<unsigned int>& tetrahedron_indices,
                  Scalar TotalMass, Scalar poisson_ratio, Scalar young_modulus)
            : bounds(generate_FEM3D_from_tetrahedron_mesh(simulation, 3 *TotalMass / tetrahedron_vertices.size(), poisson_ratio, young_modulus, tetrahedron_indices, tetrahedron_vertices)),
              TotalMass(TotalMass)
        {}

        inline Vec3 compute_center_of_mass(const Simulation& simulation, const PhysicsState& state) const {
            Vec3 center_of_mass = Vec3::Zero();
            Scalar total_mass = 0;
            for (unsigned int i = bounds.particle_index; i < bounds.particle_index + bounds.n_particles; i++) {
                const Particle& p = simulation.simulables.particles[i];
                center_of_mass += p.get_position(state.x) * p.mass;
                total_mass += p.mass;
            }
            return center_of_mass / total_mass;
        }

        inline unsigned int get_n_particles() const { return bounds.n_particles; }

        inline void freeze_particles(Simulation& simulation, const std::vector<unsigned int>& particle_indices) const {
            for (unsigned int i = 0; i < particle_indices.size(); i++) {
                const Particle& p = simulation.simulables.particles[i];
                simulation.frozen_dof.push_back(p.index+0);
                simulation.frozen_dof.push_back(p.index+1);
                simulation.frozen_dof.push_back(p.index+2);
            }
        }

    const Scalar TotalMass;
    const SimulableBounds bounds;
};

#endif // MANDOS_H_
