#ifndef FEM_UNIT_H_
#define FEM_UNIT_H_

#include "linear_algebra.hpp"
#include "particle.hpp"
#include "physics_state.hpp"

/**
 * Compute the tetrahedron linear shape function derivative with respect to the tetrahedron vertices postions.
 *
 * As the shape function is linear with the tetrahedron vertices, the derivative will be a constant matrix.
 * This derivative contains the information about the rest configuration for our FEM element.
 *
 * @param x1,x2,x3,x4 The 4 tetrahedron vertices positions.
 */
Eigen::Matrix<Scalar,4,3> compute_shape_function_derivative(const Vec3& x1, const Vec3& x2, const Vec3& x3, const Vec3& x4); // ds_dx

/**
 * Compute the deformation tensor (F) of a tetrahedron.
 *
 * @param dvecF_dx Vectorized shape function derivative. dvecF_dx = block_matrix<3,4>(ds_dx.transpose())
 * @param x1,x2,x3,x4 The 4 tetrahedron vertices positions.
 */
Mat3 compute_deformation_tensor(const Eigen::Matrix<Scalar,9,12>& dvecF_dx, const Vec3& x1, const Vec3& x2, const Vec3& x3, const Vec3& x4); // F


// FEM Materials
// Materials are structures that can compute their energy, energy gradient and hessian.
// The values are actually densities, which means that it is necessary to scale them by the
// finite element's volume.
// The energy derivatives are computed with respect to the vectorized deformation tensor vecF.
// To compute the derivatives with respect to the vertices, one must use the jacobian dvecF_dx.
// ----------------------------------------------------------------------------------------

/**
 * FEM Material
 *
 * Virtual class exposing the mandatory API that each material must implement.
 * Note that this virtual class is never instantiated at runtime, and only exists
 * for nicer compilation time errors.
 */
struct FEM_Material {
    virtual Scalar get_phi(const Mat3& F) const = 0;
    virtual Vec9 get_phi_gradient(const Mat3& F) const = 0;
    virtual Mat9 get_phi_hessian(const Mat3& F) const = 0;
};

/**
 * FEM Linear Material
 *
 * Simple material with a linear strain tensor defined as epsilon = F + F^T - I.
 * This linear strain understands rotation as deformation. The material breaks when
 * it is rotated.
 *
 * @param mu, lambda Material's Lamé coefficients
 */
struct FEM_LinearMaterial : FEM_Material {
    FEM_LinearMaterial(Scalar mu, Scalar lambda);

    // Lame coefficients
    const Scalar mu, lambda;

    Scalar get_phi(const Mat3& F) const;
    Vec9 get_phi_gradient(const Mat3& F) const;
    Mat9 get_phi_hessian(const Mat3& F) const;
};

/**
 * FEM Neo Hookean Material
 *
 * Stable Neo Hookean Material as defined in this paper https://graphics.pixar.com/library/StableElasticity/paper.pdf.
 * This material has a non linear strain which does not measure rotation as deformations.
 *
 * @param mu, lambda Material's Lamé coefficients
 */
struct FEM_NeoHookeanMaterial : FEM_Material {
    FEM_NeoHookeanMaterial(Scalar mu, Scalar lambda);

    // Lame coefficients
    const Scalar mu, lambda;

    Scalar get_phi(const Mat3& F) const;
    Vec9 get_phi_gradient(const Mat3& F) const;
    Mat9 get_phi_hessian(const Mat3& F) const;
};
// ----------------------------------------------------------------------------------------

template <typename MaterialType>
struct FEM_Element3D : PotentialEnergy {
    FEM_Element3D(Particle p1,Particle p2, Particle p3, Particle p4, Eigen::Matrix<Scalar, 4, 3> ds_dx, MaterialType material);

    const Particle p1, p2, p3, p4;
    const MaterialType material;

    // Shape function derivative
    const Eigen::Matrix<Scalar,9,12> dvecF_dx;

    Scalar compute_energy(Scalar TimeStep, const PhysicsState& state) const;

    void compute_energy_gradient(Scalar TimeStep, const PhysicsState& state, Vec& grad) const;

    void compute_energy_and_derivatives(Scalar TimeStep, const PhysicsState& state, EnergyAndDerivatives& out) const;

    Scalar compute_volume(const Vec3& x1, const Vec3& x2, const Vec3& x3, const Vec3& x4) const;
};

bool is_tetrahedron_inverted(const Vec3& v1, const Vec3& v2, const Vec3& v3, const Vec3& v4);

#endif // FEM_UNIT_H_
