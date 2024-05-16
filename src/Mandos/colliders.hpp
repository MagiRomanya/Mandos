#ifndef MANDOS_COLLIDERS_H_
#define MANDOS_COLLIDERS_H_

#include <vector>

#include <Mandos/linear_algebra.hpp>
#include <Mandos/mesh.hpp>
#include <Mandos/physics_state.hpp>

namespace tmd {
    class TriangleMeshDistance;
}

namespace mandos
{
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

// struct SDFCollider {
//     SDFCollider();
//     SDFCollider(const SDFCollider& other); // copy constructor
//     SDFCollider(const SimulationMesh& mesh);
//     ~SDFCollider();

//     SimulationMesh mesh;
//     std::unique_ptr<tmd::TriangleMeshDistance> sdf;

//     void compute_contact_geometry(const Vec3& point, ContactEvent& out) const;

// };

struct CapsuleCollider {
    Vec3 pointA, pointB;
    Scalar radius;

    void compute_contact_geometry(const Vec3& point, ContactEvent& out);
};

struct Colliders {
    std ::vector<SphereCollider> sphere_colliders;
    std ::vector<PlaneCollider> plane_colliders;
    // std ::vector<SDFCollider> sdf_colliders;

    template <typename Visitor>
    inline void for_each(Visitor visitor) const
    {
        visitor(sphere_colliders);
        visitor(plane_colliders);
        // visitor(sdf_colliders);
    }
};

struct Simulables;
void find_point_particle_contact_events(const Colliders& colliders,
                                        const Simulables& simulables,
                                        const PhysicsState& state,
                                        std::vector<ContactEvent>& events);

void compute_contact_events_energy_and_derivatives(const Scalar TimeStep,
                                                   const std::vector<ContactEvent>& events,
                                                   const PhysicsState state,
                                                   EnergyAndDerivatives& out);
}  // namespace mandos

#endif  // MANDOS_COLLIDERS_H_
