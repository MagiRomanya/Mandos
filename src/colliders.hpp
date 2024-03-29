#ifndef COLLIDERS_H_
#define COLLIDERS_H_

#include <vector>
#include "linear_algebra.hpp"
#include "physics_state.hpp"

struct ContactEvent {
    Vec3 normal;
    Scalar s_distance;
    unsigned int index;
};

struct SphereCollider {
    Vec3 center;
    Scalar radius;

    void compute_contact_geometry(const Vec3& point, ContactEvent& out) const;
};

struct PlaneCollider {
    Vec3 center, normal;

    void compute_contact_geometry(const Vec3& point, ContactEvent& out) const;
};

// struct CapsuleCollider {
//     Vec3 pointA, pointB;
//     Scalar radius;

//     void compute_contact_geometry(const Vec3& point, ContactEvent& out);
// };

#define COLLIDER_MEMBERS \
    COL(std::vector<SphereCollider>, sphere_colliders) \
    COL(std::vector<PlaneCollider>, plane_colliders)
    // COL(std::vector<CapsuleCollider>, caplule_colliders) \

#define COL(type, name) type name;

struct Colliders {
    COLLIDER_MEMBERS
};

#undef COL

struct Simulables;
void find_point_particle_contact_events(const Colliders& colliders, const Simulables& simulables, const PhysicsState& state, std::vector<ContactEvent>& events);

void compute_contact_events_energy_and_derivatives(const Scalar TimeStep, const std::vector<ContactEvent>& events, const PhysicsState state, EnergyAndDerivatives& out);

#endif // COLLIDEEnergyAndDerivativesEner
