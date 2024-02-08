#ifndef MANDOS_H_
#define MANDOS_H_

#include <cassert>
#include <vector>

#include "particle.hpp"
#include "simulable_generator.hpp"
#include "simulation.hpp"

/**
 * Rigid Body wrapper with utility functions.
 */
class RigidBodyHandle {
    public:
        RigidBodyHandle(Simulation& simulation, Scalar mass, const std::vector<Scalar> vertices);

        RigidBodyHandle(Simulation& simulation, Scalar mass, const Mat3& inertia_tensor);

        /**
         * Set the Rigid Body center of mass position in world coordinates.
         *
         * @param pos The new position of the center of mass
         */
        RigidBodyHandle set_COM_initial_position(Vec3 pos) const;

        /**
         * Set the Rigid Body initial orientation.
         *
         * @param axis_angle The new orientation using axis angle meaning angle * axis (being axis a unit vector)
         */
        RigidBodyHandle set_initial_orientation(Vec3 axis_angle) const;

        /**
         * Set the Rigid Body center of mass initial velocity.
         *
         * @param vel The new velocity of the center of mass
         */
        RigidBodyHandle set_COM_initial_velocity(Vec3 vel) const;

        /**
         * Set the Rigid Body inital angular velocity
         *
         * @param omega The angular velocity given in axis angle notation
         */
        RigidBodyHandle set_initial_angular_velocity(Vec3 omega) const;

        /**
         * Make this Rigid Body affected by gravity
         *
         * @param gravity The intensity of the gravity in the y direction
         */
        RigidBodyHandle add_gravity(Scalar gravity) const;

        /**
         * Freeze the Rigid Body degrees of freedom related to translation
         *
         */
        RigidBodyHandle freeze_translation() const;

        /**
         * Freeze the Rigid Body degrees of freedom related to rotation
         *
         * @param simulation The working simulation
         */
        RigidBodyHandle freeze_rotation() const;

        /**
         * Computes the global transformation matrix for the Rigid Body
         *
         * This method is intended to be used when rendering the rigid body.
         * Applying this transform to a certain Mesh will translate it and rotate it to
         * the proper Rigid Body position and orientation.
         *
         * @param state The global physics state in a certain instant
         */
        Mat4 get_transformation_matrix(const PhysicsState& state) const;

        const RigidBody rb;
        const unsigned int rb_index;
    private:
        Simulation& simulation;
};

class MassSpringHandle {
    public:
        MassSpringHandle(Simulation& simulation,
                         const std::vector<Scalar>& vertices,
                         const std::vector<unsigned int>& indices,
                         Scalar TotalMass, Scalar k_tension, Scalar k_bending, Scalar damping);

        Vec3 compute_center_of_mass(const PhysicsState& state) const;

        MassSpringHandle set_initial_COM_position(const Vec3& position) const;

        MassSpringHandle freeze_particles(const std::vector<unsigned int>& particle_indices) const;

        MassSpringHandle add_gravity(const Scalar gravity) const;

        inline unsigned int get_n_particles() const { return bounds.n_particles; }

        void get_dof_vector(const PhysicsState& state, std::vector<float>& out_dofs) const;

        const Scalar TotalMass;
        const SimulableBounds bounds;
    private:
        Simulation& simulation;
};

class ParticleHandle {
        public:
        ParticleHandle(Simulation& sim, const Particle& particle, unsigned int index)
                : particle(particle), particle_index(index), simulation(sim)
                {}
        ParticleHandle(Simulation& simulation, Scalar mass);

        ParticleHandle set_initial_position(Vec3 position) const;

        ParticleHandle set_initial_velocity(Vec3 velocity) const;

        ParticleHandle add_gravity(Scalar gravity) const;

        ParticleHandle freeze() const;

        inline Vec3 get_position(const PhysicsState& state) { return particle.get_position(state); }

        const Particle particle;
        const unsigned int particle_index;
    private:
        Simulation& simulation;

};

ParticleHandle get_particle_handle(Simulation& sim, unsigned int particle_index);

void join_particles_with_spring(Simulation& simulation, const ParticleHandle& p1, const ParticleHandle& p2, Scalar k, Scalar damping);

void join_rigid_body_with_particle(Simulation& sim, RigidBodyHandle rbA, ParticleHandle p);

class FEMHandle {
    public:
        FEMHandle(Simulation& simulation,
                  const std::vector<Scalar>& tetrahedron_vertices,
                  const std::vector<unsigned int>& tetrahedron_indices,
                  Scalar TotalMass, Scalar poisson_ratio, Scalar young_modulus);

        Vec3 compute_center_of_mass(const PhysicsState& state) const;

        FEMHandle freeze_particles(const std::vector<unsigned int>& particle_indices) const;

        FEMHandle add_gravity(Scalar gravity) const;

        inline unsigned int get_n_particles() const { return bounds.n_particles; }

        const Scalar TotalMass;
        const SimulableBounds bounds;
    private:
        Simulation& simulation;
};

#endif // MANDOS_H_
