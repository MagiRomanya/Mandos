#ifndef LINEAR_ALGEBRA_H_
#define LINEAR_ALGEBRA_H_

#include "../external/eigen/Eigen/Core"
#include "../external/eigen/Eigen/SparseCore"

typedef double Scalar;

// STATIC VECTORS AND MATRICES
typedef Eigen::Vector2<Scalar> Vec2;
typedef Eigen::Vector3<Scalar> Vec3;
typedef Eigen::Matrix3<Scalar> Mat3;
typedef Eigen::Vector4<Scalar> Vec4;
typedef Eigen::Matrix4<Scalar> Mat4;
typedef Eigen::Matrix4<Scalar> Mat4;
typedef Eigen::Matrix<Scalar,9,9> Mat9;
typedef Eigen::Vector<Scalar,9> Vec9;

// DYNAMIC VECTORS AND MATRICES
typedef Eigen::VectorX<Scalar> Vec;
typedef Eigen::MatrixX<Scalar> Mat;

// SPARSE MATRICES
typedef Eigen::Triplet<Scalar> Triplet;
typedef Eigen::SparseMatrix<Scalar> SparseMat;


// Given a square DxD matrix mat, return a vector of D² components
// 1  3
// 2  4  ->  vec = 1, 2, 3, 4
template<int D>
inline Eigen::Vector<Scalar, D*D> vectorize_matrix(const Eigen::Matrix<Scalar,D,D>& mat) {
  Eigen::Vector<Scalar, D*D> result = Eigen::Vector<Scalar, D*D>::Zero();
  for (unsigned int i = 0; i< D; i++) {
    for (unsigned int j = 0; j< D; j++) {
      result(D*j+i) = mat(i,j);
    }
  }
  return result;
}

template<int D>
inline Eigen::Vector<Scalar, D*D> vectorize_matrix_rowwise(const Eigen::Matrix<Scalar,D,D>& mat) {
  Eigen::Vector<Scalar, D*D> result = Eigen::Vector<Scalar, D*D>::Zero();
  for (unsigned int i = 0; i< D; i++) {
    for (unsigned int j = 0; j< D; j++) {
      result(D*j+i) = mat(j,i);
    }
  }
  return result;
}

// Given a NxM matrix, return a 3Nx3M matrix where each matrix component m_i_j becomes now m_i_j*I
// where I is the 3x3 identity matrix
// A  B            A* I, B*I
// C  D  ->  mat = C* I, D*I
template<int N, int M>
inline Eigen::Matrix<Scalar, 3*N, 3*M> block_matrix(const Eigen::Matrix<Scalar,N,M>& mat) {
  Eigen::Matrix<Scalar, 3*N, 3*M> result = Eigen::Matrix<Scalar, 3*N, 3*M>::Zero();
  for (unsigned int i = 0; i< N; i++) {
    for (unsigned int j = 0; j< M; j++) {
      // Wierd syntax we have to use because of templates
      // It just means result.block...
      result.template block<3,3>(i*3, j*3) = mat(i, j) * Mat3::Identity();
    }
  }
  return result;
}

/**
 * https://en.wikipedia.org/wiki/Levi-Civita_symbol#Three_dimensions
 */
inline Eigen::Matrix<Scalar, 3, 9> vectorized_levi_civita() {
    Eigen::Matrix<Scalar, 3, 9> e;
    e << 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f,-1.0f, 0.0f,
         0.0f, 0.0f,-1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,
         0.0f, 1.0f, 0.0f,-1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f;
    return e;
}

template<int N>
inline Eigen::Vector<Scalar, N> transpose_vectorized_vector(const Eigen::Vector<Scalar,N>& in) {
  const int n = sqrt(N);
  Eigen::Vector<Scalar, N> v;
  for (int i=0; i < n; i++) {
    for (int j=0; j < n; j++) {
      int index = n*i + j;
      int index_t = n*j + i;
      v(index) = in(index_t);
    }
  }
  return v;
}

template <int N, int M>
inline Eigen::Matrix<Scalar, N, M> transpose_vectorized_matrix(const Eigen::Matrix<Scalar, N, M>& in) {
  Eigen::Matrix<Scalar, N, M> m;
  const int n = sqrt(N);
  for (int i = 0; i < n; i++) {
    m.col(i) = transpose_vectorized_vector<N>(in.col(i));
  }
  return m;
}


#endif // LINEAR_ALGEBRA_H_
