#include <Eigen/Geometry>
#include "../simulable_generator.hpp"
#include "../rod_segment.hpp"

SimulableBounds generate_rod(Simulation& simulation,
                             unsigned int segments, Scalar TotalMass,
                             Scalar length, const Vec3 origin, const Vec3 direction,
                             const RodSegmentParameters& rod_parameters)
{
    if (segments == 0) std::cerr << "ERROR::GENERATE_ROD:: Zero segments given!" << std::endl;

    const unsigned int index = simulation.initial_state.get_nDoF();
    const unsigned int nDoF = (segments  + 1) * 6;
    const Scalar rb_mass = TotalMass / (segments + 1);

    // Resize degrees of freedom
    simulation.initial_state.add_size(nDoF);

    // Generate the rigid bodies
    // ---------------------------------------------------------------------------------
    const unsigned int rigid_body_index = simulation.simulables.rigid_bodies.size();
    const Mat3& rod_inertia_tensor = rb_mass * 10.0 * Mat3::Identity();
    for (unsigned int i = 0; i < nDoF; i+=6) {
        RigidBody rb = RigidBody(index + i, rb_mass, rod_inertia_tensor);
        add_rigid_body_to_simulation(simulation, rb);
    }

    // Rod Initial conditions
    // ---------------------------------------------------------------------------------
    const Scalar L0 = length / segments;
    const Vec3 normalized_direction = direction.normalized();

    // Compute initial rotation
    Vec3 axis_angle = Vec3::Zero();
    if (not normalized_direction.isApprox(Vec3(0.0, 0.0, 1.0), 1e-6)) {
        const Vec3 tangent = cross(normalized_direction, Vec3(0.0, 1.0, 0.0)).normalized();
        const Vec3 bitangent = cross(normalized_direction, tangent).normalized();
        Mat3 rotation;
        rotation << tangent, bitangent, normalized_direction;
        Eigen::AngleAxis<Scalar> angle_axis = Eigen::AngleAxis<Scalar>(rotation);
        axis_angle = angle_axis.axis() * angle_axis.angle();
    }

    // Finally actually setting the initial conditions.
    for (unsigned int i = 0; i < nDoF; i+=6) {
        const Vec3 position = origin + L0 * normalized_direction * (i/6);
        simulation.initial_state.x.segment<3>(index + i) = position;
        simulation.initial_state.x.segment<3>(index + i + 3) = axis_angle;
        simulation.initial_state.v.segment<6>(index + i) = Vec6::Zero();

        // // Add gravity
        // unsigned int y_index = index + i + 1;
        // simulation.energies.gravities.emplace_back(y_index, GravityParameters(-rb_mass));
    }

    // Set up the Rod Segment energies
    // ---------------------------------------------------------------------------------
    for (unsigned int i = 0; i < segments; i++) {
        RodSegmentParameters param = rod_parameters;
        param.L0 = L0;
        const RigidBody& rbA = simulation.simulables.rigid_bodies[rigid_body_index + i];
        const RigidBody& rbB = simulation.simulables.rigid_bodies[rigid_body_index + i + 1];

        simulation.energies.rod_segments.emplace_back(rbB, rbA, param);
    }

    // Compute bounds
    // ---------------------------------------------------------------------------------
    SimulableBounds bounds;
    bounds.particle_index = simulation.simulables.particles.size();
    bounds.n_particles = 0;
    bounds.nDoF = nDoF;
    bounds.n_rb = segments + 1;
    bounds.rb_index = rigid_body_index;

    return bounds;
}
