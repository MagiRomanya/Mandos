#ifndef MANDOS_LINEAR_ALGEBRA_H_
#define MANDOS_LINEAR_ALGEBRA_H_

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace mandos
{

using Scalar = double;

// STATIC VECTORS AND MATRICES
using Vec2 = Eigen::Vector2<Scalar>;
using Vec3 = Eigen::Vector3<Scalar>;
using Mat3 = Eigen::Matrix3<Scalar>;
using Vec4 = Eigen::Vector4<Scalar>;
using Mat4 = Eigen::Matrix4<Scalar>;
using Mat6 = Eigen::Matrix<Scalar, 6, 6>;
using Mat9 = Eigen::Matrix<Scalar, 9, 9>;
using Vec6 = Eigen::Vector<Scalar, 6>;
using Vec9 = Eigen::Vector<Scalar, 9>;

// DYNAMIC VECTORS AND MATRICES
using Vec = Eigen::VectorX<Scalar>;
using Mat = Eigen::MatrixX<Scalar>;

// SPARSE MATRICES
using Triplet = Eigen::Triplet<Scalar>;
using SparseMat = Eigen::SparseMatrix<Scalar, Eigen::RowMajor>;

// Given a square DxD matrix mat, return a vector of D² components
// 1  3
// 2  4  ->  vec = 1, 2, 3, 4
template <int D>
Eigen::Vector<Scalar, D * D> vectorize_matrix(const Eigen::Matrix<Scalar, D, D>& mat)
{
    // We vectorize by columns, so we can use reshaped
    return mat.reshaped();
}

template <int D>
Eigen::Vector<Scalar, D * D> vectorize_matrix_rowwise(const Eigen::Matrix<Scalar, D, D>& mat)
{
    Eigen::Vector<Scalar, D * D> result = Eigen::Vector<Scalar, D * D>::Zero();
    for (unsigned int i = 0; i < D; i++) {
        for (unsigned int j = 0; j < D; j++) {
            result(D * j + i) = mat(j, i);
        }
    }
    return result;
}

// Given a NxM matrix, return a 3Nx3M matrix where each matrix component m_i_j becomes now m_i_j*I
// where I is the 3x3 identity matrix
// A  B            A* I, B*I
// C  D  ->  mat = C* I, D*I
template <int N, int M>
Eigen::Matrix<Scalar, 3 * N, 3 * M> block_matrix(const Eigen::Matrix<Scalar, N, M>& mat)
{
    Eigen::Matrix<Scalar, 3 * N, 3 * M> result = Eigen::Matrix<Scalar, 3 * N, 3 * M>::Zero();
    for (unsigned int i = 0; i < N; i++) {
        for (unsigned int j = 0; j < M; j++) {
            // Wierd syntax we have to use because of templates
            // It just means result.block...
            result.template block<3, 3>(i * 3, j * 3) = mat(i, j) * Mat3::Identity();
        }
    }
    return result;
}

/**
 * https://en.wikipedia.org/wiki/Levi-Civita_symbol#Three_dimensions
 */
Eigen::Matrix<Scalar, 3, 9> vectorized_levi_civita();

template <int N>
Eigen::Vector<Scalar, N> transpose_vectorized_vector(const Eigen::Vector<Scalar, N>& in)
{
    const int n = static_cast<int>(sqrt(static_cast<float>(N)));
    Eigen::Vector<Scalar, N> v;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int index = n * i + j;
            int index_t = n * j + i;
            v(index) = in(index_t);
        }
    }
    return v;
}

template <int N, int M>
Eigen::Matrix<Scalar, N, M> transpose_vectorized_matrix_N(const Eigen::Matrix<Scalar, N, M>& in)
{
    Eigen::Matrix<Scalar, N, M> mat;
    const int n = static_cast<int>(sqrt(N));
    for (int i = 0; i < n; i++) {
        mat.col(i) = transpose_vectorized_vector<N>(in.col(i));
    }
    return mat;
}

template <int N, int M>
Eigen::Matrix<Scalar, N, M> transpose_vectorized_matrix_M(const Eigen::Matrix<Scalar, N, M>& in)
{
    Eigen::Matrix<Scalar, N, M> mat;
    const int m = static_cast<int>(sqrt(M));
    for (int i = 0; i < m; i++) {
        mat.row(i) = transpose_vectorized_vector<M>(in.row(i));
    }
    return mat;
}

}  // namespace mandos

#endif // MANDOS_LINEAR_ALGEBRA_H_
