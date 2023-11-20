#ifndef UTILITY_FUNCTIONS_H_
#define UTILITY_FUNCTIONS_H_

#include <raylib.h>
#include <string>
#include <vector>

#include "linear_algebra.hpp"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

Mesh LoadMeshTinyOBJ(std::string inputfile);

// Geometry related functions
Mat3 skew(const Vec3& v);
inline Scalar compute_trinagle_area(const Vec3& AB, const Vec3& AC) { return (skew(AB) * AC).norm() / 2; }
Scalar compute_tetrahedron_volume(const Vec3& AB, const Vec3& AC, const Vec3& AD);
Scalar compute_mesh_volume(const std::vector<unsigned int>& indices, const std::vector<Scalar>& vertices);
Scalar compute_mesh_surface_area(const std::vector<unsigned int>& indices, const std::vector<Scalar>& vertices);
Mat3 compute_rotation_matrix_rodrigues(const Vec3& theta);
#endif // UTILITY_FUNCTIONS_H_
