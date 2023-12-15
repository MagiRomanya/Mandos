#ifndef MANDOS_H_
#define MANDOS_H_

#include <cassert>
#include <vector>

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
         * @param simulation The working simulation
         * @param pos The new position of the center of mass
         */
        RigidBodyHandle set_COM_initial_position(Simulation& simulation, Vec3 pos) const;

        /**
         * Set the Rigid Body initial orientation.
         *
         * @param simulation The working simulation
         * @param axis_angle The new orientation using axis angle meaning angle * axis (being axis a unit vector)
         */
        RigidBodyHandle set_initial_orientation(Simulation& simulation, Vec3 axis_angle) const;

        /**
         * Set the Rigid Body center of mass initial velocity.
         *
         * @param simulation The working simulation
         * @param vel The new velocity of the center of mass
         */
        RigidBodyHandle set_COM_initial_velocity(Simulation& simulation, Vec3 vel) const;

        /**
         * Set the Rigid Body inital angular velocity
         *
         * @param simulation The working simulation
         * @param omega The angular velocity given in axis angle notation
         */
        RigidBodyHandle set_initial_angluar_velocity(Simulation& simulation, Vec3 omega) const;

        /**
         * Make this Rigid Body affected by gravity
         *
         * @param simulation The working simulation
         * @param gravity The intensity of the gravity in the y direction
         */
        RigidBodyHandle add_gravity(Simulation& simulation, Scalar gravity) const;

        /**
         * Freeze the Rigid Body degrees of freedom related to translation
         *
         * @param simulation The working simulation
         */
        RigidBodyHandle freeze_translation(Simulation& simulation) const;

        /**
         * Freeze the Rigid Body degrees of freedom related to rotation
         *
         * @param simulation The working simulation
         */
        RigidBodyHandle freeze_rotation(Simulation& simulation) const;

        /**
         * Computes the global transformation matrix for the Rigid Body
         *
         * This method is intended to be used when rendering the rigid body.
         * Applying this transform to a certain Mesh will translate it and rotate it to
         * the proper Rigid Body position and orientation.
         *
         * @param state The global physics state in a certain instant
         */
        Mat4 get_RB_transformation(const PhysicsState& state) const;

        const RigidBody rb;
        const unsigned int rb_index;
};

#ifdef ENABLE_LAGRANGE_MULTIPLIER_CONSTRAINTS
void join_rigid_bodies(Simulation& simulation, RigidBodyHandle rbA, Vec3 pA, RigidBodyHandle rbB, Vec3 pB);
#endif //ENABLE_LAGRANGE_MULTIPLIER_CONSTRAINTS


class MassSpringHandle {
    public:
        MassSpringHandle(Simulation& simulation,
                         const std::vector<Scalar>& vertices,
                         const std::vector<unsigned int>& indices,
                         Scalar TotalMass, Scalar k_tension, Scalar k_bending, Scalar damping);

        Vec3 compute_center_of_mass(const Simulation& simulation, const PhysicsState& state) const;

        MassSpringHandle set_initial_COM_position(Simulation& simulation, const Vec3& position) const;

        MassSpringHandle freeze_particles(Simulation& simulation, const std::vector<unsigned int>& particle_indices) const;

        inline unsigned int get_n_particles() const { return bounds.n_particles; }


        void get_dof_vector(const PhysicsState& state, std::vector<float>& out_dofs) const;

        const Scalar TotalMass;
    private:
        const SimulableBounds bounds;
};

class ParticleHandle {
    public:
        ParticleHandle(Simulation& simulation, Scalar mass);

        ParticleHandle set_initial_position(Simulation& simulation, Vec3 position) const;

        ParticleHandle set_initial_velocity(Simulation& simulation, Vec3 velocity) const;

        inline Vec3 get_position(const PhysicsState& state) { return particle.get_position(state.x); }

        const Particle particle;
        const unsigned int particle_index;
};

class FEMHandle {
    public:
        FEMHandle(Simulation& simulation,
                  const std::vector<Scalar>& tetrahedron_vertices,
                  const std::vector<unsigned int>& tetrahedron_indices,
                  Scalar TotalMass, Scalar poisson_ratio, Scalar young_modulus);

        Vec3 compute_center_of_mass(const Simulation& simulation, const PhysicsState& state) const;

        FEMHandle freeze_particles(Simulation& simulation, const std::vector<unsigned int>& particle_indices) const;

        inline unsigned int get_n_particles() const { return bounds.n_particles; }


        const Scalar TotalMass;
        const SimulableBounds bounds;
};

#endif // MANDOS_H_
