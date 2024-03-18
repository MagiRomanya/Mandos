#include "particle.hpp"
#include "rigid_body.hpp"


struct RodSegmentParameters {
    Scalar Ks, L0, translational_damping, rotational_damping, constraint_stiffness;
    Vec3 intrinsic_darboux, stiffness_tensor;

    Scalar get_energy(const Vec3& x1, const Vec3& x2, const Vec3& v1, const Vec3& v2,
                      const Mat3& R1, const Mat3& R2, const Mat3& R_dot1, const Mat3& R_dot2) const;

    Vec6 get_energy_gradient(Scalar TimeStep, const Vec3& x1, const Vec3& x2, const Vec3& v1, const Vec3& v2,
                             const Mat3& R1, const Mat3& R2, const Mat3& R_dot1, const Mat3& R_dot2) const;

    Mat6 get_energy_hessian(Scalar TimeStep, const Vec3& x1, const Vec3& x2, const Vec3& v1, const Vec3& v2,
                            const Mat3& R1, const Mat3& R2, const Mat3& R_dot1, const Mat3& R_dot2) const;
};

struct RodSegment {
    const Particle p1, p2;
    const RigidBody rb1, rb2;
    const RodSegmentParameters parameters;

    Scalar compute_energy(Scalar TimeStep, const PhysicsState& state) const;
    void compute_energy_and_derivatives(Scalar TimeStep, const PhysicsState& state, EnergyAndDerivatives& out) const;
};
