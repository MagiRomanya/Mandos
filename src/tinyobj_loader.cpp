#include <cassert>
#include <cstddef>
#include <limits>
#include <string>
#include <unordered_set>
#include <vector>
#include <iostream>
#include <raylib.h>

#include "linear_algebra.hpp"
#include "rigid_body.hpp"
#include "utility_functions.hpp"
#include "mesh.hpp"

#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* .cc
// Optional. define TINYOBJLOADER_USE_MAPBOX_EARCUT gives robust trinagulation. Requires C++11
// #define TINYOBJLOADER_USE_MAPBOX_EARCUT
#include "tiny_obj_loader.h"

SimulationMesh::SimulationMesh(std::string filename) {
    LoadVerticesAndIndicesTinyOBJ(filename, vertices, indices);
    recenter_mesh(*this, compute_COM_position_PARTICLES(vertices));
}

struct h_vec3 {
    h_vec3() {}
    h_vec3(Scalar x, Scalar y, Scalar z) : x(x), y(y), z(z) {}
    Scalar x, y, z;

    size_t hash() const {
        // Large prime numbers
        const size_t a = 73856093;
        const size_t b = 19349663;
        const size_t c = 83492791;
        const size_t n = std::numeric_limits<size_t>::max();
        return (((size_t)x)*a xor ((size_t)y)*b xor ((size_t)z)*c) % n;
    }

    bool operator==(const h_vec3& other) const {
        const float epsilon = 1e-6;
        return abs(x - other.x) < epsilon && abs(y - other.y) < epsilon && abs(z - other.z) < epsilon;
    }
};

namespace std {
    template <>
    struct hash<h_vec3> {
        size_t operator()(const h_vec3& v) const {
            return v.hash();
        }
    };
}

/**
 * Initialize the simulation mesh from a render mesh.
 *
 * The render mesh will have repeated vertices to properly define normals and texture
 * coordinates throughout the surface. For simulation repeated vertices are not desired
 * and are eliminated in this function.
 *
 * @param render_mesh A render mesh (only necessary to have the vertices and indices vectors initialized)
 */
SimulationMesh::SimulationMesh(const RenderMesh& render_mesh) {
    std::unordered_map<h_vec3, unsigned int> vertex_index_map; // Maps the vertex to the new simulation mesh index

    for (unsigned int i = 0; i < render_mesh.vertices.size()/9; i++) {
        // Get the triangle vertices
        const unsigned int a = 9*i + 0;
        const unsigned int b = 9*i + 3;
        const unsigned int c = 9*i + 6;
        const h_vec3 va = h_vec3(render_mesh.vertices[a], render_mesh.vertices[a+1], render_mesh.vertices[a+2]);
        const h_vec3 vb = h_vec3(render_mesh.vertices[b], render_mesh.vertices[b+1], render_mesh.vertices[b+2]);
        const h_vec3 vc = h_vec3(render_mesh.vertices[c], render_mesh.vertices[c+1], render_mesh.vertices[c+2]);

        // The new triangle indices
        unsigned int a_new;
        unsigned int b_new;
        unsigned int c_new;

        // Save only the vertices not already found
        // and initialize the new indices accordingly
        if (not vertex_index_map.contains(va)) {
            a_new = vertices.size()/3;
            vertex_index_map[va] = a_new;
            vertices.push_back(va.x);
            vertices.push_back(va.y);
            vertices.push_back(va.z);
        }
        else {
            a_new = vertex_index_map[va];
        }
        if (not vertex_index_map.contains(vb)) {
            b_new = vertices.size()/3;
            vertex_index_map[vb] = b_new;
            vertices.push_back(vb.x);
            vertices.push_back(vb.y);
            vertices.push_back(vb.z);
        }
        else {
            b_new = vertex_index_map[vb];
        }
        if (not vertex_index_map.contains(vc)) {
            c_new = vertices.size()/3;
            vertex_index_map[vc] = c_new;
            vertices.push_back(vc.x);
            vertices.push_back(vc.y);
            vertices.push_back(vc.z);
        }
        else {
            c_new = vertex_index_map[vc];
        }

        // Store the new triangle indices
        indices.push_back(a_new);
        indices.push_back(b_new);
        indices.push_back(c_new);
    }
}

void LoadRenderMeshTinyOBJ(std::string inputfile,
                                   std::vector<float> &out_vertices,
                                   std::vector<float> &out_normals,
                                   std::vector<float> &out_tex_coord) {
    tinyobj::ObjReaderConfig reader_config;
    tinyobj::ObjReader reader;
    if (!reader.ParseFromFile(inputfile, reader_config)) {
        if (!reader.Error().empty()) {
            std::cerr << "ERROR::TinyObjReader: " << reader.Error();
        } else {
            std::cerr << "ERROR::TinyObjReader: unknown cause :v" << std::endl;
        }
        exit(1);
    }

    if (!reader.Warning().empty()) {
        std::cout << "WARNING::TinyObjReader: " << reader.Warning();
    }

    auto &attrib = reader.GetAttrib();
    auto &shapes = reader.GetShapes();

    // Loop over shapes
    for (size_t s = 0; s < shapes.size(); s++) {
        // Loop over faces(polygon)
        size_t index_offset = 0;
        for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
            size_t fv = size_t(shapes[s].mesh.num_face_vertices[f]);

            // Loop over vertices in the face.
            for (size_t v = 0; v < fv; v++) {
              // access to vertex
              tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
              tinyobj::real_t vx =
                  attrib.vertices[3 * size_t(idx.vertex_index) + 0];
              tinyobj::real_t vy =
                  attrib.vertices[3 * size_t(idx.vertex_index) + 1];
              tinyobj::real_t vz =
                  attrib.vertices[3 * size_t(idx.vertex_index) + 2];
              out_vertices.push_back(vx);
              out_vertices.push_back(vy);
              out_vertices.push_back(vz);

              // Check if `normal_index` is zero or positive. negative = no normal data
              if (idx.normal_index >= 0) {
                tinyobj::real_t nx = attrib.normals[3 * size_t(idx.normal_index) + 0];
                tinyobj::real_t ny = attrib.normals[3 * size_t(idx.normal_index) + 1];
                tinyobj::real_t nz = attrib.normals[3 * size_t(idx.normal_index) + 2];
                out_normals.push_back(nx);
                out_normals.push_back(ny);
                out_normals.push_back(nz);
              }

              // Check if `texcoord_index` is zero or positive. negative = no
              // texcoord data
              if (idx.texcoord_index >= 0) {
                tinyobj::real_t tx =
                    attrib.texcoords[2 * size_t(idx.texcoord_index) + 0];
                tinyobj::real_t ty =
                    attrib.texcoords[2 * size_t(idx.texcoord_index) + 1];
                out_tex_coord.push_back(tx);
                    out_tex_coord.push_back(ty);
              }
            }
            index_offset += fv;
        }
    }
}

inline Vec3 compute_tangent_vector(const Vec3& edge1, const Vec3& edge2, const Vec2& deltaUV1, const Vec2& deltaUV2) {
    Vec3 tangent;
    const float denominator = 1.0f / (deltaUV1.x() * deltaUV2.y() - deltaUV1.y() * deltaUV2.x());
    // First vertex tangent
    tangent.x() = denominator * (deltaUV2.y() * edge1.x() - deltaUV1.y() * edge2.x());
    tangent.y() = denominator * (deltaUV2.y() * edge1.y() - deltaUV1.y() * edge2.y());
    tangent.z() = denominator * (deltaUV2.y() * edge1.z() - deltaUV1.y() * edge2.z());
    return tangent.normalized();
}

inline void computeRenderMeshTangentVectors(RenderMesh& mesh) {
    mesh.tangents.resize(mesh.vertices.size());
    const unsigned int nTriangles = mesh.vertices.size() / 9;
    for (unsigned int i = 0; i < nTriangles; ++i) {
        // Triangle edges
        const Vec3 x1 = Vec3(mesh.vertices[9*i + 3*0 + 0], mesh.vertices[9*i + 3*0 + 1], mesh.vertices[9*i + 3*0 + 2]);
        const Vec3 x2 = Vec3(mesh.vertices[9*i + 3*1 + 0], mesh.vertices[9*i + 3*1 + 1], mesh.vertices[9*i + 3*1 + 2]);
        const Vec3 x3 = Vec3(mesh.vertices[9*i + 3*2 + 0], mesh.vertices[9*i + 3*2 + 1], mesh.vertices[9*i + 3*2 + 2]);

        // Texture coordinates
        const Vec2 uv1 = Vec2(mesh.texcoords[6*i + 2*0 + 0], mesh.texcoords[6*i + 2*0 + 1]);
        const Vec2 uv2 = Vec2(mesh.texcoords[6*i + 2*1 + 0], mesh.texcoords[6*i + 2*1 + 1]);
        const Vec2 uv3 = Vec2(mesh.texcoords[6*i + 2*2 + 0], mesh.texcoords[6*i + 2*2 + 1]);

        const Vec3 edge1 = x2 - x1;
        const Vec3 edge2 = x3 - x1;
        const Vec2 deltaUV1 = uv2 - uv1;
        const Vec2 deltaUV2 = uv3 - uv1;
        const Vec3 tangent = compute_tangent_vector(edge1, edge2, deltaUV1, deltaUV2);

        // First vertex tangent
        mesh.tangents[9 * i + 3 * 0 + 0] = tangent.x();
        mesh.tangents[9 * i + 3 * 0 + 1] = tangent.y();
        mesh.tangents[9 * i + 3 * 0 + 2] = tangent.z();

        // Second vertex tangent
        mesh.tangents[9 * i + 3 * 1 + 0] = tangent.x();
        mesh.tangents[9 * i + 3 * 1 + 1] = tangent.y();
        mesh.tangents[9 * i + 3 * 1 + 2] = tangent.z();

        // Third vertex tangent
        mesh.tangents[9 * i + 3 * 2 + 0] = tangent.x();
        mesh.tangents[9 * i + 3 * 2 + 1] = tangent.y();
        mesh.tangents[9 * i + 3 * 2 + 2] = tangent.z();
    }
}

RenderMesh::RenderMesh(std::string filename) {
    LoadRenderMeshTinyOBJ(filename, vertices, normals, texcoords);
    computeRenderMeshTangentVectors(*this);
}

void LoadVerticesAndIndicesTinyOBJ(std::string inputfile, std::vector<float>& out_vertices, std::vector<unsigned int>& out_indices) {
    tinyobj::ObjReaderConfig reader_config;
    tinyobj::ObjReader reader;
    if (!reader.ParseFromFile(inputfile, reader_config)) {
        if (!reader.Error().empty()) {
            std::cerr << "TinyObjReader: " << reader.Error();
        }
        else {
            std::cerr << "TinyObjReader: unknown cause"<< std::endl;
        }
        exit(1);
    }

    auto& attrib = reader.GetAttrib();
    auto& shapes = reader.GetShapes();

    out_vertices = attrib.vertices;

    for (size_t s = 0; s < shapes.size(); s++) {
        const tinyobj::shape_t& shape = shapes[s];
        for (size_t i = 0; i < shape.mesh.indices.size(); i++){
            size_t index = shape.mesh.indices[i].vertex_index;
            out_indices.push_back(index);
        }
    }
}

void recenter_mesh(SimulationMesh& mesh, const Vec3& com) {
    for (unsigned int i = 0; i < mesh.vertices.size() / 3; i++) {
        mesh.vertices[3*i+0] -= com.x();
        mesh.vertices[3*i+1] -= com.y();
        mesh.vertices[3*i+2] -= com.z();
    }
}

void RenderMesh::updateFromSimulationMesh(const SimulationMesh& sim_mesh) {
    // The number of triangles of both meshes must be the same
    assert(sim_mesh.indices.size() / 3 == mesh.vertices.size() / 9);
    if (tangents.size() == 0)
        tangents.resize(vertices.size());

    const unsigned int nTriangles = sim_mesh.indices.size() / 3;
    for (unsigned int i = 0; i < nTriangles; i++) {
        // each triangle has 3 mesh.vertices -> 9 coordinates
        const Vec3 x1 = Vec3(sim_mesh.vertices[3*sim_mesh.indices[3*i+0] + 0], sim_mesh.vertices[3*sim_mesh.indices[3*i+0] + 1] , sim_mesh.vertices[3*sim_mesh.indices[3*i+0] + 2]);
        const Vec3 x2 = Vec3(sim_mesh.vertices[3*sim_mesh.indices[3*i+1] + 0], sim_mesh.vertices[3*sim_mesh.indices[3*i+1] + 1] , sim_mesh.vertices[3*sim_mesh.indices[3*i+1] + 2]);
        const Vec3 x3 = Vec3(sim_mesh.vertices[3*sim_mesh.indices[3*i+2] + 0], sim_mesh.vertices[3*sim_mesh.indices[3*i+2] + 1] , sim_mesh.vertices[3*sim_mesh.indices[3*i+2] + 2]);

        // Compute the normal vector ( the same for all the mesh.vertices in the triangle )
        const Vec3 edge1 = (x2 - x1).normalized();
        const Vec3 edge2 = (x3 - x1).normalized();
        const Vec3 normal = cross(edge1, edge2);

        // Compute the tangent vector ( the same for all the mesh.vertices in the triangle )
        const Vec2 uv1 = Vec2(texcoords[6*i + 2*0 + 0], texcoords[6*i + 2*0 + 1]);
        const Vec2 uv2 = Vec2(texcoords[6*i + 2*1 + 0], texcoords[6*i + 2*1 + 1]);
        const Vec2 uv3 = Vec2(texcoords[6*i + 2*2 + 0], texcoords[6*i + 2*2 + 1]);
        const Vec2 deltaUV1 = uv2 - uv1;
        const Vec2 deltaUV2 = uv3 - uv1;
        const Vec3 tangent = compute_tangent_vector(edge1, edge2, deltaUV1, deltaUV2);

        for (unsigned int j = 0; j < 3; j++) {
            const unsigned int index = 3*sim_mesh.indices[3*i+j];
            // Update the mesh.vertices
            vertices[9*i + 3*j + 0] = sim_mesh.vertices[index + 0]; // x
            vertices[9*i + 3*j + 1] = sim_mesh.vertices[index + 1]; // y
            vertices[9*i + 3*j + 2] = sim_mesh.vertices[index + 2]; // z

            // Update the normals
            normals[9*i + 3*j + 0] = normal.x();
            normals[9*i + 3*j + 1] = normal.y();
            normals[9*i + 3*j + 2] = normal.z();

            // Update the tangents
            tangents[9*i + 3*j + 0] = tangent.x();
            tangents[9*i + 3*j + 1] = tangent.y();
            tangents[9*i + 3*j + 2] = tangent.z();
        }

    }

}
