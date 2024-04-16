#include <Eigen/Geometry>
#include "../simulable_generator.hpp"
#include "../rod_segment.hpp"

Vec3 compute_axis_angle_from_direction(const Vec3& direction) {
    const Vec3 normalized_direction = direction.normalized();

    // Compute initial rotation
    if (not normalized_direction.isApprox(Vec3(0.0, 0.0, 1.0), 1e-6)) {
        const Vec3 tangent = cross(normalized_direction, Vec3(0.0, 1.0, 0.0)).normalized();
        const Vec3 bitangent = cross(normalized_direction, tangent).normalized();
        Mat3 rotation;
        rotation << tangent, bitangent, normalized_direction;
        Eigen::AngleAxis<Scalar> angle_axis = Eigen::AngleAxis<Scalar>(rotation);
        return angle_axis.axis() * angle_axis.angle();
    }
    return Vec3(0.0, 0.0, 0.0);
}


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
    const unsigned int dof_index = static_cast<unsigned int>(simulation.initial_state.x.size());
    simulation.initial_state.add_size(nDoF);

    // Generate the rigid bodies
    // ---------------------------------------------------------------------------------
    const unsigned int rigid_body_index = static_cast<unsigned int>(simulation.simulables.rigid_bodies.size());
    const Mat3& rod_inertia_tensor = rb_mass * 10.0 * Vec3(1.0, 1.0, 2.0).asDiagonal();
    for (unsigned int i = 0; i < nDoF; i+=6) {
        RigidBody rb = RigidBody(index + i, rb_mass, rod_inertia_tensor);
        add_rigid_body_to_simulation(simulation, rb);
    }

    // Rod Initial conditions
    // ---------------------------------------------------------------------------------
    const Scalar L0 = length / segments;
    const Vec3 normalized_direction = direction.normalized();
    const Vec3 axis_angle = compute_axis_angle_from_direction(direction);

    // Finally actually setting the initial conditions.
    for (unsigned int i = 0; i < nDoF; i+=6) {
        const Vec3 position = origin + L0 * normalized_direction * (i/6);
        simulation.initial_state.x.segment<3>(index + i) = position;
        simulation.initial_state.x.segment<3>(index + i + 3) = axis_angle;
        simulation.initial_state.v.segment<6>(index + i) = Vec6::Zero();
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
    bounds.particle_index = static_cast<unsigned int>(simulation.simulables.particles.size());
    bounds.n_particles = 0;
    bounds.nDoF = nDoF;
    bounds.n_rb = segments + 1;
    bounds.rb_index = rigid_body_index;
    bounds.dof_index = dof_index;

    return bounds;
}

SimulableBounds generate_rod(Simulation& simulation,
                             std::vector<Scalar> vertices,
                             Scalar TotalMass,
                             const RodSegmentParameters& rod_parameters)

{

    if (vertices.size() < 6) std::cerr << "ERROR::GENERATE_ROD:: Not enough vertices!" << std::endl;
    if (std::fmod(vertices.size(), 3) != 0) std::cerr << "ERROR::GENERATE_ROD:: Vertices not formated correctly!" << std::endl;

    const unsigned int index = simulation.initial_state.get_nDoF();
    const unsigned int nRB = static_cast<unsigned int>(vertices.size()) / 3;
    const unsigned int nSegments = nRB - 1;
    const unsigned int nDoF =  nRB * 6;
    const Scalar rb_mass = TotalMass / nRB;

    // Resize degrees of freedom
    const unsigned int dof_index = static_cast<unsigned int>(simulation.initial_state.x.size());
    simulation.initial_state.add_size(nDoF);

    // Generate the rigid bodies
    // ---------------------------------------------------------------------------------
    const unsigned int rigid_body_index = static_cast<unsigned int>(simulation.simulables.rigid_bodies.size());
    const Mat3& rod_inertia_tensor = rb_mass * 10.0 * Vec3(1.0, 1.0, 2.0).asDiagonal();
    for (unsigned int i = 0; i < nDoF; i+=6) {
        RigidBody rb = RigidBody(index + i, rb_mass, rod_inertia_tensor);
        add_rigid_body_to_simulation(simulation, rb);
    }

    // Rod Initial conditions
    // ---------------------------------------------------------------------------------
    // Compute directions
    // The first and last rigid bodies
    const Vec3 v0 = Vec3(vertices[0], vertices[1], vertices[2]);
    const Vec3 v1 = Vec3(vertices[3], vertices[4], vertices[5]);
    const Vec3 direction0 = v1 - v0;
    const Vec3 axis_angle0 = compute_axis_angle_from_direction(direction0);
    simulation.initial_state.x.segment<3>(index) = v0;
    simulation.initial_state.x.segment<3>(index+3) = axis_angle0;
    simulation.initial_state.v.segment<3>(index) = Vec3::Zero();
    simulation.initial_state.v.segment<3>(index+3) = Vec3::Zero();

    const unsigned int last_idx = static_cast<unsigned int>(vertices.size()) - 6;
    const Vec3 vLL = Vec3(vertices[last_idx], vertices[last_idx + 1], vertices[last_idx + 2]);
    const Vec3 vL = Vec3(vertices[last_idx + 3], vertices[last_idx + 4], vertices[last_idx + 5]);
    const Vec3 directionL = vL - vLL;
    const Vec3 axis_angleL = compute_axis_angle_from_direction(directionL);

    simulation.initial_state.x.segment<3>(index + nDoF - 6) = vL;
    simulation.initial_state.x.segment<3>(index + nDoF - 6 + 3) = axis_angleL;
    simulation.initial_state.v.segment<3>(index + nDoF - 6) = Vec3::Zero();
    simulation.initial_state.v.segment<3>(index + nDoF - 6 + 3) = Vec3::Zero();

    for (unsigned int i = 1; i < nSegments; i++) {
        const Vec3 vA = Vec3(vertices[3*(i-1)+0], vertices[3*(i-1)+1], vertices[3*(i-1)+2]);
        const Vec3 vB = Vec3(vertices[3*(i+0)+0], vertices[3*(i+0)+1], vertices[3*(i+0)+2]);
        const Vec3 vC = Vec3(vertices[3*(i+1)+0], vertices[3*(i+1)+1], vertices[3*(i+1)+2]);

        const Vec3 AB = (vB - vA);
        const Vec3 BC = (vC - vB);
        const Vec3 direction = 0.5 * (AB + BC);
        const Vec3 axis_angle = compute_axis_angle_from_direction(direction);

        simulation.initial_state.x.segment<3>(index + 6 * i) = vB;
        simulation.initial_state.x.segment<3>(index + 6 * i + 3) = axis_angle;
        simulation.initial_state.v.segment<6>(index + 6 * i) = Vec6::Zero();
    }

    // Set up the Rod Segment energies
    // ---------------------------------------------------------------------------------
    for (unsigned int i = 0; i < nSegments; i++) {
        const RigidBody& rbA = simulation.simulables.rigid_bodies[rigid_body_index + i];
        const RigidBody& rbB = simulation.simulables.rigid_bodies[rigid_body_index + i + 1];
        RodSegmentParameters param = rod_parameters;
        const PhysicsState& i_state = simulation.initial_state;
        param.L0 = (rbA.get_COM_position(i_state.x) - rbB.get_COM_position(i_state.x)).norm();
        param.intrinsic_darboux = compute_darboux_vector(param.L0, rbB.compute_rotation_matrix(i_state.x), rbA.compute_rotation_matrix(i_state.x));
        simulation.energies.rod_segments.emplace_back(rbB, rbA, param);
    }

    // Compute bounds
    // ---------------------------------------------------------------------------------
    SimulableBounds bounds;
    bounds.particle_index = static_cast<unsigned int>(simulation.simulables.particles.size());
    bounds.n_particles = 0;
    bounds.nDoF = nDoF;
    bounds.n_rb = nRB;
    bounds.rb_index = rigid_body_index;
    bounds.dof_index = dof_index;

    return bounds;
}
