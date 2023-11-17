#include <vector>

#include "physics_state.hpp"
#include "rigid_body.hpp"
#include "linear_algebra.hpp"
#include "utility_functions.hpp"

Mat3 skew(const Vec3& v) {
    Mat3 m;
    m <<     0.,  -v.z(),   v.y(),
          v.z(),      0.,  -v.x(),
         -v.y(),   v.x(),      0.;
    return m;
}

Mat3 compute_rotation_matrix_rodrigues(const Vec3& theta) {
    const Scalar angle = theta.norm();

    // Prevent angle = 0, axis not well defined case
    // Using small angle approximation instead
    const Scalar THRESHOLD = 1e-1;
    if (angle < THRESHOLD) return Mat3::Identity() + skew(theta);

    const Vec3 axis = theta / angle;
    // Rodrigues formula https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula#Matrix_notation
    const Mat3 K = skew(axis);
    return Mat3::Identity() + sin(angle) * K + (1. -cos(angle)) * K * K;
}

Mat3 RigidBody::compute_inertia_tensor(const Mat3& rotation_matrix) const {
    return rotation_matrix * inertia_tensor0 * rotation_matrix.transpose();
}

Scalar RigidBody::get_kinetic_energy(const Vec3& v, const Vec3& omega, const Mat3& inertia_tensor) const {
    const Scalar rectiliniar_energy = 0.5 * mass * v.squaredNorm();
    const Scalar rotational_energy = 0.5 * omega.transpose() * inertia_tensor * omega;
    return rectiliniar_energy + rotational_energy;
}

Vec3 RigidBody::get_coriolis_torque(const Vec3& omega, const Mat3& inertia_tensor) const {
    const Vec3 coriolis_torque = - skew(omega) * inertia_tensor * omega;
    return coriolis_torque;
}

void RigidBody::compute_energy_and_derivatives(const PhysicsState& state, EnergyAndDerivatives& out) const {
    // Get the relevant sate
    // ---------------------------------------------------------------
    const Vec3 x = state.x.segment(index, 3);
    const Vec3 v = state.v.segment(index, 3);
    const Vec3 theta = state.x.segment(index + 3, 3);
    const Vec3 omega = state.v.segment(index + 3, 3);

    // Compute rotation matrix and Inerta tensor
    // ---------------------------------------------------------------
    const Mat3 rotation_matrix = compute_rotation_matrix_rodrigues(theta);
    const Mat3 inertia_tensor = compute_inertia_tensor(rotation_matrix);

    // Compute the energy derivatives
    // ---------------------------------------------------------------
    const Scalar kinetic_energy = get_kinetic_energy(v, omega, inertia_tensor);
    const Vec3 coriolis_torque = get_coriolis_torque(omega, inertia_tensor);

    // Add the energy derivatives to the global structure
    // ---------------------------------------------------------------

    out.energy += kinetic_energy;
    for (unsigned int i = 0; i<3; i++) {
        out.force[index + 3 + i] += coriolis_torque(i);  // torque
        // TODO: add force derivatives
        // for (unsigned int j = 0; j<3; j++) {
        //     out.df_dx_triplets.push_back(Triplet(index+i, index+j, df_dx(i, j)));
        // }
    }
}

Vec3 RigidBody::get_COM_position(const PhysicsState& state) const {
    return state.x.segment(index, 3);
}

Mat3 RigidBody::compute_rotation_matrix(const PhysicsState& state) const {
    Vec3 theta = state.x.segment(index+3, 3);
    return compute_rotation_matrix_rodrigues(theta);
}

Scalar compute_tetrahedron_volume(const Vec3& AB, const Vec3& AC, const Vec3& AD) {
    return (skew(AB) * AC).dot(AD) / 6.;
}

// Computethe volume of the mesh in the same units as the vertex positions.
// This function assumes the mesh to be convex
Scalar compute_mesh_volume(const std::vector<unsigned int>& indices, const std::vector<Scalar>& vertices) {
    Vec3 mesh_center = Vec3::Zero();
    Scalar volume = 0;
    // Compute mesh center
    // ---------------------------------------------------------------
    for (unsigned int v=0; v < vertices.size(); v+=3) {
        mesh_center += Vec3(vertices[v], vertices[v+1], vertices[v+2]);
    }
    mesh_center /= vertices.size() / 3.0;

    // Iterate over triangles and add the volume of each tetrahedron
    // ---------------------------------------------------------------
    for (unsigned int t=0; t < indices.size(); t+=3) {
        const Vec3 v1 = Vec3(vertices[3*indices[t]], vertices[3*indices[t]+1], vertices[3*indices[t]+2]);
        const Vec3 v2 = Vec3(vertices[3*indices[t+1]], vertices[3*indices[t+1]+1], vertices[3*indices[t+1]+2]);
        const Vec3 v3 = Vec3(vertices[3*indices[t+2]], vertices[3*indices[t+2]+1], vertices[3*indices[t+2]+2]);
        volume += abs(compute_tetrahedron_volume(v1 - mesh_center, v2 - mesh_center, v3 - mesh_center));
    }
    return volume;
}

Scalar compute_mesh_surface_area(const std::vector<unsigned int>& indices, const std::vector<Scalar>& vertices) {
    Scalar surface_area = 0.0;
    for (unsigned int t=0; t < indices.size(); t+=3) {
        const Vec3 v1 = Vec3(vertices[3*indices[t]], vertices[3*indices[t]+1], vertices[3*indices[t]+2]);
        const Vec3 v2 = Vec3(vertices[3*indices[t+1]], vertices[3*indices[t+1]+1], vertices[3*indices[t+1]+2]);
        const Vec3 v3 = Vec3(vertices[3*indices[t+2]], vertices[3*indices[t+2]+1], vertices[3*indices[t+2]+2]);
        surface_area += abs(compute_trinagle_area(v2 - v1, v3 - v1));
    }
    return surface_area;
}

Mat3 compute_initial_inertia_tensor(Scalar rb_total_mass, const std::vector<unsigned int>& indices, const std::vector<Scalar>& vertices, RB_MASS_DISTRIBUTION mass_distribution) {
    Mat3 inertia_tensor = Mat3::Zero();
    switch (mass_distribution) {
        case PARTICLES:
            {
                const Scalar PARTICLE_MASS = rb_total_mass * 3.0 / vertices.size();
                for (unsigned int i = 0; i < vertices.size(); i+=3) {
                    const Vec3 r = Vec3(vertices[i], vertices[i+1], vertices[i+2]);
                    const Mat3 r_star = skew(r);
                    inertia_tensor += - PARTICLE_MASS * r_star * r_star;
                }
            }
            break;
        case SHELL:
            {

            }
            break;
        case UNIFORM_VOLUME:
            {

            }
            break;
    }
    return inertia_tensor;
}

Vec3 compute_COM_position(const std::vector<unsigned int>& indices, const std::vector<Scalar>& vertices, RB_MASS_DISTRIBUTION mass_distribution) {
    Vec3 COM_position = Vec3::Zero();
    switch (mass_distribution) {
        case PARTICLES:
        {
            for (unsigned int i = 0; i < vertices.size(); i+=3) {
                const Vec3 r = Vec3(vertices[i], vertices[i+1], vertices[i+2]);
                COM_position +=r;
            }
            // Average the result
            COM_position *= 3.0 / vertices.size();
        }
        break;
        case SHELL:
        {
            Scalar total_surface = 0.0;
            for (unsigned int t=0; t < indices.size(); t+=3) {
                const Vec3 v1 = Vec3(vertices[3*indices[t]], vertices[3*indices[t]+1], vertices[3*indices[t]+2]);
                const Vec3 v2 = Vec3(vertices[3*indices[t+1]], vertices[3*indices[t+1]+1], vertices[3*indices[t+1]+2]);
                const Vec3 v3 = Vec3(vertices[3*indices[t+2]], vertices[3*indices[t+2]+1], vertices[3*indices[t+2]+2]);
                const Scalar triangle_area = abs(compute_trinagle_area(v2 - v1, v3 - v1));
                COM_position += triangle_area * (v1 + v2 + v3) / 3;
                total_surface += total_surface;
            }
            // Average the result
            COM_position /= total_surface;
        }
        case UNIFORM_VOLUME:
        {
            const Vec3 origin = compute_COM_position(indices, vertices, PARTICLES);
            Scalar total_volume = 0.0;
            for (unsigned int t=0; t < indices.size(); t+=3) {
                const Vec3 v1 = Vec3(vertices[3*indices[t]], vertices[3*indices[t]+1], vertices[3*indices[t]+2]);
                const Vec3 v2 = Vec3(vertices[3*indices[t+1]], vertices[3*indices[t+1]+1], vertices[3*indices[t+1]+2]);
                const Vec3 v3 = Vec3(vertices[3*indices[t+2]], vertices[3*indices[t+2]+1], vertices[3*indices[t+2]+2]);
                const Scalar tetrahedron_volume = abs(compute_tetrahedron_volume(v1-origin, v2-origin, v3-origin));
                total_volume += tetrahedron_volume;
                COM_position += tetrahedron_volume * (v1 + v2 + v3 + origin) / 4;
            }
            COM_position /= total_volume;
        }
    }
    return COM_position;
}

Mat3 RigidBody::compute_inertia_tensor(const PhysicsState& state) const {
    return compute_inertia_tensor(compute_rotation_matrix(state));
}
