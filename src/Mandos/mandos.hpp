#ifndef MANDOS_H_
#define MANDOS_H_

#include <cassert>
#include <vector>

#include <Mandos/simulable_generator.hpp>
#include <Mandos/simulation.hpp>

/**
 * Here there are the definitions of simulable handles, that are a front end of the simulation for the user to interact with.
 *
 * The handle objects proide the user with a friendly way of creating simulables and creating interactions and constraints between them.
 * This abstracts the user from the implementation details of the creation process of simulabes, that have to allocate degrees of freedom and construct the simulables from atomic energies.
 *
 * The main idea here is that in order to create the simulable, the user must pass a simulation instance where it wants the simulable to be created.
 * All the important state is thus in the simulation instance and the simulable Handle only serves as an easy way to edit the simulation instance.
 * NO REAL STATE IS STORED IN THE HANDLES.
 */

/**
 * Rigid Body wrapper with utility functions.
 */
class RigidBodyHandle {
    public:
        RigidBodyHandle(Simulation& simulation, Scalar mass, const std::vector<Scalar> vertices, bool global = false);

        RigidBodyHandle(Simulation& simulation, Scalar mass, const Mat3& inertia_tensor, bool global = false);

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

        Scalar compute_energy(const PhysicsState& state, const PhysicsState& state0);

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

void join_rigid_body_com_with_spring(Simulation& simulation, const RigidBodyHandle& rbA, const RigidBodyHandle& rbB, Scalar k, Scalar damping);

void join_rigid_body_with_spring(Simulation& simulation, const RigidBodyHandle& rbA, const Vec3& pA, const RigidBodyHandle& rbB, const Vec3& pB, Scalar k, Scalar damping);

void join_rigid_body_with_rod_segment(Simulation& simulation, RigidBodyHandle& rbA, RigidBodyHandle& rbB, RodSegmentParameters parameters);

class FEMHandle {
    public:
        FEMHandle(Simulation& simulation,
                  const std::vector<Scalar>& tetrahedron_vertices,
                  const std::vector<unsigned int>& tetrahedron_indices,
                  Scalar TotalMass, Scalar poisson_ratio, Scalar young_modulus);

        Vec3 compute_center_of_mass(const PhysicsState& state) const;

        FEMHandle freeze_particles(const std::vector<unsigned int>& particle_indices) const;

        FEMHandle add_gravity(Scalar gravity) const;

        FEMHandle set_com_position(Vec3 position) const;

        inline unsigned int get_n_particles() const { return bounds.n_particles; }

        const Scalar TotalMass;
        const SimulableBounds bounds;
    private:
        Simulation& simulation;
};

class RodHandle {
    public:
        RodHandle(Simulation& simulation, unsigned int segments, const Vec3& rod_vector, Scalar TotalMass, const RodSegmentParameters& parameters);

        RodHandle(Simulation& simulation, const std::vector<Scalar>& vertices, Scalar TotalMass, const RodSegmentParameters& parameters);

        Vec3 compute_center_of_mass(const PhysicsState& state) const;

        RodHandle add_gravity(Scalar gravity) const;
        RodHandle set_rigid_body_initial_position(const Scalar s, const Vec3& x) const;
        RodHandle set_rigid_body_initial_velocity(const Scalar s, const Vec3& v) const;
        RodHandle set_rigid_body_initial_orientation(const Scalar s, const Vec3& x) const;
        RodHandle set_rigid_body_initial_angular_velocity(const Scalar s, const Vec3& v) const;
        RodHandle set_initial_rod_position(const Vec3& origin) const;
        RodHandle freeze_rigid_body(Scalar s) const;

        inline unsigned int get_n_rigid_bodies() const { return bounds.n_rb; }

        const Scalar TotalMass;
        const SimulableBounds bounds;
    private:
        Simulation& simulation;
};

class PlaneColliderHandle {
    public:
        PlaneColliderHandle(Simulation& simulation);

        PlaneColliderHandle set_origin_position(const Vec3& origin) const;
        PlaneColliderHandle set_direction(const Vec3& direction) const;

        unsigned int index;

    private:
        Simulation& simulation;
};

class SphereColliderHandle {
    public:
        SphereColliderHandle(Simulation& simulation);

        SphereColliderHandle set_origin_position(const Vec3& origin) const;
        SphereColliderHandle set_radius(const Scalar& radius) const;

        unsigned int index;

    private:
        Simulation& simulation;
};

class SDFColliderHandle {
    public:
        SDFColliderHandle(Simulation& simulation, const SimulationMesh& sim_mesh);

        unsigned int index;
    private:
        Simulation& simulation;
};

#endif // MANDOS_H_
