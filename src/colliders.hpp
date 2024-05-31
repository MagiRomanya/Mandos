#ifndef COLLIDERS_H_
#define COLLIDERS_H_

#include <memory>
#include <vector>
#include "linear_algebra.hpp"
#include "mesh.hpp"
#include "physics_state.hpp"
#include "utility_functions.hpp"

namespace tmd {
    class TriangleMeshDistance;
}

struct ContactEvent {
    Vec3 normal;
    Scalar s_distance;
    unsigned int index;
};

struct Collider {
    virtual void compute_contact_geometry(const Vec3& point, ContactEvent& out) const = 0;
};

struct SphereCollider final : Collider {
    Vec3 center;
    Scalar radius;

    void compute_contact_geometry(const Vec3& point, ContactEvent& out) const;
};

struct PlaneCollider final : Collider {
    Vec3 center, normal;

    void compute_contact_geometry(const Vec3& point, ContactEvent& out) const;
};

struct SDFCollider final : Collider {
    SDFCollider();
    SDFCollider(const SDFCollider& other); // copy constructor
    SDFCollider(const SimulationMesh& mesh);
    ~SDFCollider();

    SimulationMesh mesh;
    std::unique_ptr<tmd::TriangleMeshDistance> sdf;

    void compute_contact_geometry(const Vec3& point, ContactEvent& out) const;
};

// struct CapsuleCollider {
//     Vec3 pointA, pointB;
//     Scalar radius;

//     void compute_contact_geometry(const Vec3& point, ContactEvent& out);
// };

struct Colliders {
    std::vector<SphereCollider> sphere_colliders;
    std::vector<PlaneCollider> plane_colliders;
    std::vector<SDFCollider> sdf_colliders;

    template <typename Visitor>
    constexpr inline void for_each(Visitor && visitor) const {
        visitor(sphere_colliders);
        visitor(plane_colliders);
        visitor(sdf_colliders);
    }
};

CHECK_WETHER_COMPOSITE_IS_VALID(Colliders, Collider);

struct Simulables;
void find_point_particle_contact_events(const Colliders& colliders, const Simulables& simulables, const PhysicsState& state, std::vector<ContactEvent>& events);

void compute_contact_events_energy_and_derivatives(const Scalar TimeStep, const std::vector<ContactEvent>& events, const PhysicsState state, EnergyAndDerivatives& out);

#endif // COLLIDERS_H_