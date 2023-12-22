#include <cstddef>
#include <limits>
#include <string>
#include <unordered_set>
#include <vector>
#include <iostream>
#include <raylib.h>

#include "linear_algebra.hpp"
#include "utility_functions.hpp"
#include "mesh.hpp"

#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* .cc
// Optional. define TINYOBJLOADER_USE_MAPBOX_EARCUT gives robust trinagulation. Requires C++11
// #define TINYOBJLOADER_USE_MAPBOX_EARCUT
#include "tiny_obj_loader.h"

Mesh LoadMeshFromVectors(std::vector<unsigned short>& indices, const std::vector<float>& vertices, const std::vector<float>& normals, const std::vector<float>& texcoords) {
    Mesh mesh = {0};
    mesh.vertexCount = vertices.size() / 3;
    mesh.triangleCount = indices.size() / 3;
    mesh.vertices = (float *)std::memcpy(RL_MALLOC(vertices.size() * sizeof(float)), vertices.data(), vertices.size() * sizeof(float));
    mesh.texcoords = (float *)std::memcpy(RL_MALLOC(texcoords.size() * sizeof(float)), texcoords.data(), texcoords.size() * sizeof(float));
    mesh.normals = (float *)std::memcpy(RL_MALLOC(normals.size() * sizeof(float)), normals.data(), normals.size() * sizeof(float));
    mesh.indices = (unsigned short *)std::memcpy(RL_MALLOC(indices.size() * sizeof(unsigned short)), indices.data(), indices.size() * sizeof(unsigned short));
    UploadMesh(&mesh, false);
    return mesh;
}

Mesh LoadMeshTinyOBJ(std::string inputfile) {
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

    std::vector<unsigned short> indices;
    std::vector<float> vertices;
    std::vector<float> normals;
    std::vector<float> texcoords;

    texcoords.resize(attrib.vertices.size() / 3 * 2, 0);
    normals.resize(attrib.vertices.size());

    for (size_t s = 0; s < shapes.size(); s++) {
        const tinyobj::shape_t& shape = shapes[s];
        for (size_t i = 0; i < shape.mesh.indices.size(); i++){
            size_t index = shape.mesh.indices[i].vertex_index;
            indices.push_back(index);

            // Align the vertex attributes with the vertices, so that vertex indeces will
            // also point to the correct vertex attributes.
            // NOTE This is not perfect as it assigns one and only one attribute to each vertex,
            // which can translate to not ideal normals and texture coordinates.
            size_t uv_index = shape.mesh.indices[i].texcoord_index;
            texcoords[index*2] = attrib.texcoords[2*uv_index];
            texcoords[index*2+1] = attrib.texcoords[2*uv_index+1];

            size_t normal_index = shape.mesh.indices[i].normal_index;
            normals[index*3] = attrib.normals[3*normal_index];
            normals[index*3+1] = attrib.normals[3*normal_index+1];
            normals[index*3+2] = attrib.normals[3*normal_index+2];
        }
    }

    vertices = attrib.vertices;
    return LoadMeshFromVectors(indices, vertices, normals, texcoords);
}


SimulationMesh::SimulationMesh(std::string filename) {
    LoadVerticesAndIndicesTinyOBJ(filename, vertices, indices);
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

    for (unsigned int i = 0; i < render_mesh.indices.size()/3; i++) {
        // Get the triangle vertices
        const unsigned int a = 3*render_mesh.indices[3*i+0];
        const unsigned int b = 3*render_mesh.indices[3*i+1];
        const unsigned int c = 3*render_mesh.indices[3*i+2];
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

void LoadVerticesAndIndicesTinyOBJ(std::string inputfile,
                                   std::vector<float> &out_vertices,
                                   std::vector<unsigned int> &out_indices,
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
    auto &materials = reader.GetMaterials();

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

              // Check if `normal_index` is zero or positive. negative = no
              // normal data
              if (idx.normal_index >= 0) {
                tinyobj::real_t nx =
                    attrib.normals[3 * size_t(idx.normal_index) + 0];
                tinyobj::real_t ny =
                    attrib.normals[3 * size_t(idx.normal_index) + 1];
                tinyobj::real_t nz =
                    attrib.normals[3 * size_t(idx.normal_index) + 2];
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
              out_indices.push_back(out_indices.size());
            }
            index_offset += fv;
        }
    }
}

RenderMesh::RenderMesh(std::string filename) {
    LoadVerticesAndIndicesTinyOBJ(filename, vertices, indices, normals, texcoord);
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

Mesh RenderMesh_to_RaylibMesh(const RenderMesh& render_mesh) {
    std::vector<unsigned short> indices = std::vector<unsigned short>(render_mesh.indices.begin(), render_mesh.indices.end());
    return LoadMeshFromVectors(indices, render_mesh.vertices, render_mesh.normals, render_mesh.texcoord);
}


Mesh SimulationMesh_to_RaylibMesh(const SimulationMesh& sim_mesh) {
    std::vector<unsigned short> indices = std::vector<unsigned short>(sim_mesh.indices.begin(), sim_mesh.indices.end());
    Mesh mesh = {0};
    mesh.vertexCount = sim_mesh.vertices.size() / 3;
    mesh.triangleCount = indices.size() / 3;
    mesh.vertices = (float *)std::memcpy(RL_MALLOC(sim_mesh.vertices.size() * sizeof(float)), sim_mesh.vertices.data(), sim_mesh.vertices.size() * sizeof(float));
    mesh.texcoords = NULL;
    mesh.normals = NULL;
    mesh.indices = (unsigned short *)std::memcpy(RL_MALLOC(indices.size() * sizeof(unsigned short)), indices.data(), indices.size() * sizeof(unsigned short));
    UploadMesh(&mesh, false);
    return mesh;
}
