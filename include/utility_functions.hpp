#ifndef UTILITY_FUNCTIONS_H_
#define UTILITY_FUNCTIONS_H_

#include <string>
#include <vector>
#include <iostream>
#include "linear_algebra.hpp"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#define DEBUG_LOG(variable) std::cout << #variable << " " << variable << std::endl;

/**
 * Compute the skew-symetric matrix of the given vector.
 */
Mat3 skew(const Vec3& v);

/**
 * Compute the area of the triangle.
 *
 * @param AB, AC the 2 vectors defining the triangle.
 */
inline Scalar compute_trinagle_area(const Vec3& AB, const Vec3& AC) { return (skew(AB) * AC).norm() / 2; }

/**
 * Compute the volume of the tetrahedron.
 *
 * @param AB, AC, AD the 3 vectors defining the tetrahedron
 */
Scalar compute_tetrahedron_volume(const Vec3& AB, const Vec3& AC, const Vec3& AD);

/**
 * Computes the rotation matrix from a given axis-angle rotation vector using Rodrigues'.
 */
Mat3 compute_rotation_matrix_rodrigues(const Vec3& theta);

inline Vec3 cross(const Vec3& v, const Vec3& u) {
    return Vec3(v.y() * u.z() - v.z() * u.y(), v.z() * u.x() - v.x() * u.z(), v.x() * u.y() - v.y() * u.x());
}
#endif // UTILITY_FUNCTIONS_H_
