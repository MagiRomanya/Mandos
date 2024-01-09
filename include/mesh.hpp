#ifndef MESH_H_
#define MESH_H_

#include <vector>
#include "edge.hpp"
#include "linear_algebra.hpp"

struct SimulationMesh;

/*
** Mesh description with repeated vertices, along with normals and texture coordinates.
**
** The vertices are stored as a vector of floats. The vector has size nVertices * 3. The number of triangles is nVertices / 3.
** The normals and texture coordinates are also stored as vector of floats with sizes nVertices*3 and nVertices*2 respectivly.
 */
struct RenderMesh {
    RenderMesh();
    RenderMesh(std::string filename);

    void updateFromSimulationMesh(const SimulationMesh& sim_mesh);

    std::vector<float> vertices;
    std::vector<float> normals;
    std::vector<float> texcoord;
};

/*
** Mesh description with no repeated vertices.
**
** The vertices are stored as a vector of floats. The vector has size nVertices * 3.
** The indices are stored as a vector of unsigned integers. The vector has size nTriangles * 3.
 */
struct SimulationMesh {
    SimulationMesh() {};
    SimulationMesh(std::string filename);
    SimulationMesh(const RenderMesh& render_mesh);

    std::vector<Scalar> vertices;
    std::vector<unsigned int> indices;
};

/**
 * Load an obj file to a vector of vertices and a vector of indices.
 *
 * @param inputfile directory of the obj file
 * @param out_vertices output vector of vertices (no repetition)
 * @param out_indices output vector of indices
*/
void LoadVerticesAndIndicesTinyOBJ(std::string inputfile, std::vector<float>& out_vertices, std::vector<unsigned int>& out_indices);

/**
 * Compute a tetrahedron volumetric mesh from a triangle surface mesh.
 *
 * @param triangle_indices, triangle_vertices description of the mesh with no repeated vertices
 * @param out_tetrahedron_indices Tetrahedron output vector of indices (mutiple of 4)
 * @param out_tetrahedron_vertices Tetrahedron output vector of vertices
 */
void tetgen_compute_tetrahedrons(const std::vector<unsigned int>& triangle_indices, const std::vector<float>& triangle_vertices,
                                 std::vector<unsigned int>& out_tetrahedron_indices, std::vector<float>& out_tetrahedron_vertices);

/**
 * Compute the volume of the mesh (units of the vertices).
 *
 * @param indices, vertices description of the mesh (no vertex repetition)
 */
Scalar compute_mesh_volume(const std::vector<unsigned int>& indices, const std::vector<Scalar>& vertices);

/**
 * Compute the surface area of the mesh (units of the vertices).
 *
 * @param indices, vertices description of the mesh (no vertex repetition)
 */
Scalar compute_mesh_surface_area(const std::vector<unsigned int>& indices, const std::vector<Scalar>& vertices);

/**
 * Compute the boundary of a given mesh.
 *
 * @param indices, vertices description of the mesh (no vertex repetition)
 * @param internalEdges output vector of edges inside of the mesh.
 * @param externalEdges output vector of edges in the boundary of the mesh. (no external edges means a closed mesh)
 */
void mesh_boundary(const std::vector<Scalar>& vertices,
                   const std::vector<unsigned int>& indices,
                   std::vector<Edge> &internalEdges,
                   std::vector<Edge> &externalEdges);


std::array<unsigned int, 2> count_springs(const std::vector<Scalar>& vertices, const std::vector<unsigned int>& indices);

void recenter_mesh(SimulationMesh& mesh, const Vec3& com);

#endif // MESH_H_
