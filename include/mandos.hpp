#ifndef MANDOS_H_
#define MANDOS_H_

#include <vector>

#include "gravity.hpp"
#include "linear_algebra.hpp"
#include "inertia_energies.hpp"
#include "physics_state.hpp"
#include "rigid_body.hpp"
#include "simulation.hpp"


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


#endif // MANDOS_H_
