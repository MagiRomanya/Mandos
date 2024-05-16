#include <Mandos/spring.hpp>
#include <Mandos/linear_algebra.hpp>
#include <Mandos/utility_functions.hpp>

namespace mandos
{

inline Scalar get_particle_spring_energy(const SpringParameters& param, Scalar L)
{
    return 0.5 * param.k * (L - param.L0) * (L - param.L0);
}

inline Vec3 get_particle_spring_energy_gradient(const SpringParameters& param,
                                                const Vec3& x1,
                                                const Vec3& x2,
                                                const Vec3& v1,
                                                const Vec3& v2,
                                                Scalar L)
{
    /* Computes the spring force */
    const Vec3 u = (x1 - x2) / L;
    Vec3 f = -param.k * (L - param.L0) * u;
    // damping force
    f += -param.damping * u * u.transpose() * (v1 - v2);

    // The gradient is minus the force
    return -f;
}

inline Mat3 get_particle_spring_energy_hessian(const SpringParameters& param,
                                               Scalar TimeStep,
                                               const Vec3& x1,
                                               const Vec3& x2,
                                               Scalar L)
{
    // u is the normalized vector between particles 1 and 2
    const Vec3 u = (x1 - x2) / L;
    const Mat3 uut = u * u.transpose();

    // Initialize the derivative matrix
    Mat3 df_dx = -param.k / L * ((L - param.L0) * Mat3::Identity() + param.L0 * uut);

    // Damping jacobian
    df_dx += -1.0 / TimeStep * param.damping * uut;

    // The hessian is minus the force jacobian
    return -df_dx;  // 3x3 matrix
}

Scalar ParticleSpring::compute_energy(Scalar TimeStep, const PhysicsState& state) const
{
    const Vec3 x1 = p1.get_position(state);
    const Vec3 x2 = p2.get_position(state);
    const Scalar L = (x1 - x2).norm();

    const Scalar energy = get_particle_spring_energy(parameters, L);
    return energy;
}

void ParticleSpring::compute_energy_gradient(Scalar TimeStep, const PhysicsState& state, Vec& grad) const
{
    // Get the relevant sate
    // ---------------------------------------------------------------
    const Vec3 x1 = p1.get_position(state);
    const Vec3 x2 = p2.get_position(state);
    const Vec3 v1 = p1.get_velocity(state);
    const Vec3 v2 = p2.get_velocity(state);
    const Scalar epsilon = 1e-8;
    const Scalar L = (x1 - x2).norm() + epsilon;  // Avoid division by zero

    const Vec3 gradient = get_particle_spring_energy_gradient(parameters, x1, x2, v1, v2, L);

    grad.segment<3>(p1.index) += gradient;
    grad.segment<3>(p2.index) += -gradient;
}

void ParticleSpring::compute_energy_and_derivatives(Scalar TimeStep,
                                                    const PhysicsState& state,
                                                    EnergyAndDerivatives& out) const
{
    // Get the relevant sate
    // ---------------------------------------------------------------
    const Vec3 x1 = p1.get_position(state);
    const Vec3 x2 = p2.get_position(state);
    const Vec3 v1 = p1.get_velocity(state);
    const Vec3 v2 = p2.get_velocity(state);
    const Scalar epsilon = 1e-8;
    const Scalar L = (x1 - x2).norm() + epsilon;  // Avoid division by zero

    // Compute the energy derivatives
    // ---------------------------------------------------------------
    const Scalar energy = get_particle_spring_energy(parameters, L);
    const Vec3 gradient = get_particle_spring_energy_gradient(parameters, x1, x2, v1, v2, L);
    const Mat3 hessian = get_particle_spring_energy_hessian(parameters, TimeStep, x1, x2, L);

    // Add the energy derivatives to the global structure
    // ---------------------------------------------------------------
    out.energy += energy;
    // Newton's third law: equal and opposite reaction
    out.gradient.segment<3>(p1.index) += gradient;
    out.gradient.segment<3>(p2.index) += -gradient;
    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            out.hessian_triplets.push_back(Triplet(p1.index + i, p1.index + j, hessian(i, j)));
            out.hessian_triplets.push_back(Triplet(p1.index + i, p2.index + j, -hessian(i, j)));
            out.hessian_triplets.push_back(Triplet(p2.index + i, p1.index + j, -hessian(i, j)));
            out.hessian_triplets.push_back(Triplet(p2.index + i, p2.index + j, hessian(i, j)));
        }
    }
}

inline Scalar get_rigid_body_spring_energy(const SpringParameters& param, const Scalar L)
{
    return 0.5 * (L - param.L0) * (L - param.L0);
}

inline Vec6 get_rigid_body_spring_gradientA(const SpringParameters& param,
                                            const Vec3& u,
                                            const Vec3& posA,
                                            const Mat3& RA,
                                            const Scalar L)
{
    Vec6 grad;

    const Vec3 base = param.k * (L - param.L0) * u;
    grad.segment<3>(0) = base;
    grad.segment<3>(3) = base.transpose() * skew(-RA * posA);

    return grad;
}

inline Mat6 get_rigid_body_spring_hessAA(const SpringParameters& param,
                                         const Vec3& u,
                                         const Vec3& posA,
                                         const Mat3& RA,
                                         const Scalar L)
{
    Mat6 hess;
    const Mat3 uut = u * u.transpose();
    // Linear block
    hess.block<3, 3>(0, 0) = param.k / L * ((L - param.L0) * Mat3::Identity() + param.L0 * uut);

    const Mat3 jac = skew(-RA * posA);
    // Off diagonal terms
    hess.block<3, 3>(0, 3) = hess.block<3, 3>(0, 0) * jac;
    hess.block<3, 3>(3, 0) = hess.block<3, 3>(0, 3).transpose();

    // Rotation block
    hess.block<3, 3>(3, 3) = param.k * jac.transpose() * uut * jac;
    // hess.block<3,3>(3, 3) = param.k * jac.transpose() * jac;

    return hess;
}

inline Mat6 get_rigid_body_spring_hessAB(const SpringParameters& param,
                                         const Vec3& u,
                                         const Vec3& posA,
                                         const Mat3& RA,
                                         const Vec3& posB,
                                         const Mat3& RB,
                                         const Scalar L)
{
    Mat6 hess;
    const Mat3 uut = u * u.transpose();
    // Linear block
    hess.block<3, 3>(0, 0) = -param.k / L * ((L - param.L0) * Mat3::Identity() + param.L0 * uut);

    const Mat3 jacA = skew(-RA * posA);
    const Mat3 jacB = skew(-RB * posB);
    // Off diagonal terms
    hess.block<3, 3>(0, 3) = hess.block<3, 3>(0, 0) * jacA;  // dTdx
    hess.block<3, 3>(3, 0) = -param.k * uut * jacA;          // dFdt

    // Rotation block
    hess.block<3, 3>(3, 3) = param.k * jacB.transpose() * uut * jacA;

    return hess;
}

Scalar RigidBodySpring::compute_energy(Scalar TimeStep, const PhysicsState& state) const
{
    // Get the relevant sate
    // ---------------------------------------------------------------
    const Vec3 comA = rbA.get_COM_position(state.x);
    const Vec3 comB = rbB.get_COM_position(state.x);
    const Mat3 RA = rbA.compute_rotation_matrix(state.x);
    const Mat3 RB = rbB.compute_rotation_matrix(state.x);

    const Vec3 worldA = RA * posA + comA;
    const Vec3 worldB = RB * posB + comB;
    const Scalar L = (worldA - worldB).norm();
    const Scalar energy = get_rigid_body_spring_energy(parameters, L);
    return energy;
}

void RigidBodySpring::compute_energy_gradient(Scalar TimeStep, const PhysicsState& state, Vec& grad) const
{
    // Get the relevant sate
    // ---------------------------------------------------------------
    const Vec3 comA = rbA.get_COM_position(state.x);
    const Vec3 comB = rbB.get_COM_position(state.x);
    const Mat3 RA = rbA.compute_rotation_matrix(state.x);
    const Mat3 RB = rbB.compute_rotation_matrix(state.x);

    // Avoid repeated computations
    // ---------------------------------------------------------------
    const Vec3 worldA = RA * posA + comA;
    const Vec3 worldB = RB * posB + comB;
    const Scalar L = (worldA - worldB).norm();
    const Vec3 u = (worldA - worldB) / L;

    const Vec6 gradA = get_rigid_body_spring_gradientA(parameters, u, posA, RA, L);
    const Vec6 gradB = get_rigid_body_spring_gradientA(parameters, -u, posB, RB, L);

    grad.segment<6>(rbA.index) += gradA;
    grad.segment<6>(rbB.index) += gradB;
}

void RigidBodySpring::compute_energy_and_derivatives(Scalar TimeStep,
                                                     const PhysicsState& state,
                                                     EnergyAndDerivatives& out) const
{
    // Get the relevant sate
    // ---------------------------------------------------------------
    const Vec3 comA = rbA.get_COM_position(state.x);
    const Vec3 comB = rbB.get_COM_position(state.x);
    const Mat3 RA = rbA.compute_rotation_matrix(state.x);
    const Mat3 RB = rbB.compute_rotation_matrix(state.x);

    // Avoid repeated computations
    // ---------------------------------------------------------------
    const Vec3 worldA = RA * posA + comA;
    const Vec3 worldB = RB * posB + comB;
    const Scalar L = (worldA - worldB).norm();
    const Vec3 u = (worldA - worldB) / L;
    // DEBUG_LOG(u.transpose());

    // Compute energy and derivatives
    // ---------------------------------------------------------------
    const Scalar energy = get_rigid_body_spring_energy(parameters, L);

    const Vec6 gradA = get_rigid_body_spring_gradientA(parameters, u, posA, RA, L);
    const Vec6 gradB = get_rigid_body_spring_gradientA(parameters, -u, posB, RB, L);

    Mat6 hessAA_f = Mat6::Zero();
    const Scalar dx = 1e-4;
    for (unsigned int i = 0; i < 3; i++) {
        Vec3 vdx = Vec3::Zero();
        vdx(i) = dx;
        const Vec3 dcomA = comA + vdx;
        const Mat3 dRA = compute_rotation_matrix_rodrigues(vdx) * RA;

        const Vec3 dworldAx = RA * posA + dcomA;
        const Scalar dLx = (dworldAx - worldB).norm();
        const Vec3 dux = (dworldAx - worldB) / dLx;

        const Vec3 dworldAR = dRA * posA + comA;
        const Scalar dLR = (dworldAR - worldB).norm();
        const Vec3 duR = (dworldAR - worldB) / dLx;

        const Vec6 dgradAx = get_rigid_body_spring_gradientA(parameters, dux, posA, RA, dLx);
        const Vec6 dgradAR = get_rigid_body_spring_gradientA(parameters, duR, posA, dRA, dLR);
        hessAA_f.col(i) = (dgradAx - gradA) / dx;
        hessAA_f.col(3 + i) = (dgradAR - gradA) / dx;
    }

    const Mat6 hessAA = get_rigid_body_spring_hessAA(parameters, u, posA, RA, L);
    const Mat6 hessAB = get_rigid_body_spring_hessAB(parameters, u, posA, RA, posB, RB, L);
    const Mat6 hessBB = get_rigid_body_spring_hessAA(parameters, -u, posB, RB, L);

    // std ::cout << "hessAA"
    //            << "\n" << hessAA << std ::endl;
    // std ::cout << "hessAA_f"
    //            << "\n" << hessAA_f << std ::endl;

    // Add the energy derivatives to the global structure
    // ---------------------------------------------------------------
    out.energy += energy;
    out.gradient.segment<6>(rbA.index) += gradA;
    out.gradient.segment<6>(rbB.index) += gradB;
    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            out.hessian_triplets.push_back(Triplet(rbA.index + i, rbA.index + j, hessAA(i, j)));
            out.hessian_triplets.push_back(Triplet(rbA.index + i, rbB.index + j, hessAB(i, j)));
            out.hessian_triplets.push_back(Triplet(rbB.index + i, rbA.index + j, hessAB(j, i)));
            out.hessian_triplets.push_back(Triplet(rbB.index + i, rbB.index + j, hessBB(i, j)));
        }
    }
}

}  // namespace mandos
