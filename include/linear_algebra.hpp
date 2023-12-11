#ifndef LINEAR_ALGEBRA_H_
#define LINEAR_ALGEBRA_H_

#include "../external/eigen/Eigen/Core"
#include "../external/eigen/Eigen/SparseCore"

typedef float Scalar;

// STATIC VECTORS AND MATRICES
typedef Eigen::Vector3<Scalar> Vec3;
typedef Eigen::Matrix3<Scalar> Mat3;
typedef Eigen::Vector4<Scalar> Vec4;
typedef Eigen::Matrix4<Scalar> Mat4;
typedef Eigen::Matrix4<Scalar> Mat4;

// DYNAMIC VECTORS AND MATRICES
typedef Eigen::VectorX<Scalar> Vec;
typedef Eigen::MatrixX<Scalar> Mat;

// SPARSE MATRICES
typedef Eigen::Triplet<Scalar> Triplet;
typedef Eigen::SparseMatrix<Scalar> SparseMat;


#endif // LINEAR_ALGEBRA_H_
