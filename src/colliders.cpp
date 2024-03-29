#include "colliders.hpp"
#include "simulation.hpp"


void SphereCollider::compute_contact_geometry(const Vec3& point, ContactEvent& out) const {
    const Vec3 delta = point - center;
    const Scalar distance = delta.norm();
    out.normal = delta / distance;
    out.s_distance = distance - radius;
}

void PlaneCollider::compute_contact_geometry(const Vec3& point, ContactEvent& out) const {
    out.s_distance = normal.dot(point - center);
    out.normal = normal;
}

void find_point_particle_contact_events(const Colliders& colliders, const Simulables& simulables, const PhysicsState& state, std::vector<ContactEvent>& events) {
#define COL(type, name)                                               \
    for (unsigned int j = 0; j < colliders.name.size(); j++) {        \
        colliders.name[j].compute_contact_geometry(point, event);     \
        event.index = index;                                          \
        if ( event.s_distance < 0.0 ) events.push_back(event);        \
    }

    for (unsigned int i = 0; i < simulables.particles.size(); i++) {
        const Vec3 point = simulables.particles[i].get_position(state);
        const unsigned int index = simulables.particles[i].index;
        ContactEvent event;
        COLLIDER_MEMBERS
    }

    for (unsigned int i = 0; i < simulables.rigid_bodies.size(); i++) {
        const Vec3 point = simulables.rigid_bodies[i].get_COM_position(state.x);
        const unsigned int index = simulables.rigid_bodies[i].index;
        ContactEvent event;
        COLLIDER_MEMBERS
    }

#undef COL
}

void compute_contact_events_energy_and_derivatives(const Scalar TimeStep, const std::vector<ContactEvent>& events, const PhysicsState state, EnergyAndDerivatives& out) {
    const Scalar contact_stiffness = 10000.0;
    const Scalar friction_coefficient = 100.0;
    for (unsigned int i = 0; i < events.size(); i++) {
        // Compute the energy and derivatives
        const ContactEvent& event = events[i];
        // Collision penalty force
        const Scalar c_energy = contact_stiffness * event.s_distance * event.s_distance;
        const Vec3 c_gradient = contact_stiffness * event.normal * event.s_distance;
        const Mat3 c_hessian = contact_stiffness * event.normal * event.normal.transpose();

        // Friction
        const Vec3 velocity = state.v.segment<3>(event.index);
        const Mat3 proj_mat = (Mat3::Identity() - event.normal * event.normal.transpose());
        const Vec3 tan_vel = proj_mat * velocity;
        const Scalar f_energy = 0.5 * friction_coefficient * tan_vel.squaredNorm();
        const Vec3 f_gradient = friction_coefficient * tan_vel;
        const Mat3 f_hessian = 1.0 / TimeStep * friction_coefficient * proj_mat;

        // Fill the global structure
        out.energy += c_energy + f_energy;

        out.gradient.segment<3>(event.index) += c_gradient + f_gradient;

        for (unsigned int a = 0; a < 3; a++)
            for (unsigned int b = 0; b < 3; b++)
                out.hessian_triplets.emplace_back(event.index + a, event.index + b, c_hessian(a,b) + f_hessian(a,b));
    }
}
