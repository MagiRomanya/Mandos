#ifndef UTILITY_FUNCTIONS_H_
#define UTILITY_FUNCTIONS_H_

#include <raylib.h>
#include <string>
#include <vector>
#include <iostream>

#include "linear_algebra.hpp"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#define DEBUG_LOG(variable) std::cout << #variable << " " << variable << std::endl;

Mesh LoadMeshTinyOBJ(std::string inputfile);

void LoadVerticesAndIndicesTinyOBJ(std::string inputfile, std::vector<float>& out_vertices, std::vector<unsigned int>& out_indices);

void tetgen_compute_tetrahedrons(const std::vector<unsigned int>& triangle_indices, const std::vector<float>& triangle_vertices,
                                 std::vector<unsigned int>& out_tetrahedron_indices, std::vector<float>& out_tetrahedron_vertices);

// Geometry related functions
Mat3 skew(const Vec3& v);

inline Scalar compute_trinagle_area(const Vec3& AB, const Vec3& AC) { return (skew(AB) * AC).norm() / 2; }

Scalar compute_tetrahedron_volume(const Vec3& AB, const Vec3& AC, const Vec3& AD);

Scalar compute_mesh_volume(const std::vector<unsigned int>& indices, const std::vector<Scalar>& vertices);

Scalar compute_mesh_surface_area(const std::vector<unsigned int>& indices, const std::vector<Scalar>& vertices);

Mat3 compute_rotation_matrix_rodrigues(const Vec3& theta);

inline Matrix matrix_eigen_to_raylib(const Mat4& m) {
    Matrix r = {
    m(0,0), m(0,1), m(0, 2), m(0, 3),
    m(1,0), m(1,1), m(1, 2), m(1, 3),
    m(2,0), m(2,1), m(2, 2), m(2, 3),
    m(3,0), m(3,1), m(3, 2), m(3, 3),
    };
    return r;
}

inline Vector3 vector3_eigen_to_raylib(const Vec3& v) {
    return Vector3{v.x(), v.y(), v.z()};
}


#endif // UTILITY_FUNCTIONS_H_
