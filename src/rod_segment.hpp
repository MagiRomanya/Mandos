#ifndef ROD_SEGMENT_H_
#define ROD_SEGMENT_H_

#include "rigid_body.hpp"


Vec3 compute_darboux_vector(const Scalar L0, const Mat3& R1, const Mat3& R2);

Mat3 compute_darboux_vector_local_derivative(const Scalar L0, const Mat3& R1, const Mat3& R2);

/**
 * Precompute these values so that the compute energy does not need to do it every time
 */
struct RodSegmentPrecomputedValues {
    RodSegmentPrecomputedValues() {}
    RodSegmentPrecomputedValues(Scalar L0, Scalar TimeStep,
                                const Vec3& x1, const Vec3& x2,
                                const Vec3& v1, const Vec3& v2,
                                const Mat3& R1, const Mat3& R2,
                                const Mat3& R_dot1, const Mat3& R_dot2);
    Scalar one_over_L0, one_over_h;
    Vec3 x1, x2, v1, v2;
    Vec3 v_rel;
    Vec3 deltaX;
    Scalar L, one_over_L;
    Vec3 darboux_vector;
    Mat3 darboux_vector_derivative;
    Vec3 u;
    Mat3 uut;
    Mat3 R1, R2, R_dot1, R_dot2, R;
    Vec3 C;
};

struct RodSegmentParameters {
    Scalar Ks, L0, translational_damping, rotational_damping, constraint_stiffness;
    Vec3 intrinsic_darboux, stiffness_tensor;

    Scalar compute_energy(const RodSegmentPrecomputedValues& values) const;

    Vec3 compute_energy_linear_gradient(const RodSegmentPrecomputedValues& values) const;

    Vec3 compute_energy_rotational_gradient_A(const RodSegmentPrecomputedValues& values) const;

    Vec3 compute_energy_rotational_gradient_B(const RodSegmentPrecomputedValues& values) const;

    Mat6 compute_energy_hessian_A(const RodSegmentPrecomputedValues& values) const;

    Mat6 compute_energy_hessian_B(const RodSegmentPrecomputedValues& values) const;

    Mat6 compute_energy_hessian_AB(const RodSegmentPrecomputedValues& values) const;
};

struct RodSegment {
    RodSegment(const RigidBody& rb1,  const RigidBody& rb2, const RodSegmentParameters& parameters)
        : rbA(rb1), rbB(rb2), parameters(parameters) {}

    const RigidBody rbA, rbB;
    const RodSegmentParameters parameters;

    Scalar compute_energy(Scalar TimeStep, const PhysicsState& state) const;
    void compute_energy_and_derivatives(Scalar TimeStep, const PhysicsState& state, EnergyAndDerivatives& out) const;
};

#endif // ROD_SEGMENT_H_
