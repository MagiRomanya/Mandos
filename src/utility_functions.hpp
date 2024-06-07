#ifndef UTILITY_FUNCTIONS_H_
#define UTILITY_FUNCTIONS_H_

#include <string>
#include <type_traits>
#include <vector>
#include <iostream>
#include "linear_algebra.hpp"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#define DEBUG_LOG(variable) std::cout << #variable << " " << variable << std::endl

/**
 * Compute the skew-symetric matrix of the given vector.
 */
Mat3 skew(const Vec3& v);

Vec3 unskew(const Mat& m);

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
 * Compute the rotation matrix from a given axis-angle rotation vector using Rodrigues'.
 *
 * @param theta Rotation vector defined as angle * axis, where angle is a scalar and axis a normalized vector.
 */
Mat3 compute_rotation_matrix_rodrigues(const Vec3& theta);

/**
 * Compute the cross product between 2 vectors.
 */
inline Vec3 cross(const Vec3& v, const Vec3& u) {
    return Vec3(v.y() * u.z() - v.z() * u.y(),
                v.z() * u.x() - v.x() * u.z(),
                v.x() * u.y() - v.y() * u.x());
}

template <typename T>
constexpr size_t sizeof_composite() {
    size_t size = 0;
    T instance;
    instance.for_each([&](const auto& members){
        size += sizeof(members);
    });
    return size;
}

template <typename T, typename P>
constexpr bool is_base_of_composite() {
    T instance;
    bool result = true;
    instance.for_each([&](auto vectors){
        using vector_value_type = typename decltype(vectors)::value_type;
        static_assert(std::is_base_of<P, vector_value_type>::value);
        if constexpr (not std::is_base_of<P, vector_value_type>::value) {
            result = false;
        }
    });
    return result;
}

#define CHECK_WETHER_COMPOSITE_IS_VALID(composite_type, abstract_type) \
    static_assert(sizeof_composite<composite_type>() == sizeof(composite_type), \
                  "The for_each method in" #composite_type " does not consider all members."); \
    static_assert(is_base_of_composite<composite_type, abstract_type>(), \
                  "Not all types in the composite " #composite_type " inherit from the abstract type " #abstract_type);

#endif // UTILITY_FUNCTIONS_H_
