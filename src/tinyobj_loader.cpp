#include <cassert>
#include <cstddef>
#include <string>
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
    indices = render_mesh.indices;
    const unsigned int nVertices = 1 + *std::max_element(indices.begin(), indices.end());
    const unsigned int nTriangles = static_cast<unsigned int>(indices.size()) / 3;
    vertices.resize(3*nVertices, 0.0);
    for (unsigned int i = 0; i < nTriangles; i++) {
        for (unsigned int j = 0; j < 3; j++) { // Iterate the 3 vertices of the triangle i
            const unsigned int idx = 3*indices[3*i+j];
            const Vec3 vert = Vec3(render_mesh.vertices[9*i + 0 + 3*j], render_mesh.vertices[9*i + 1 + 3*j], render_mesh.vertices[9*i + 2 + 3*j]);
            // Here we overwrite the vertices that are repeated in the Render Mesh,
            // resulting in a simulation mesh with no repeated vertices!
            vertices[idx + 0] = vert.x();
            vertices[idx + 1] = vert.y();
            vertices[idx + 2] = vert.z();
        }
    }
    recenter_mesh(*this, compute_COM_position_PARTICLES(vertices));
}

void LoadRenderMeshTinyOBJ(std::string inputfile,
                           std::vector<Scalar> &out_vertices,
                           std::vector<Scalar> &out_normals,
                           std::vector<Scalar> &out_tex_coord,
                           std::vector<unsigned int> &out_indices) {
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
        std::cout << "WARNING::TinyObjReader: " << reader.Warning() << std::endl;
    }

    auto &attrib = reader.GetAttrib();
    auto &shapes = reader.GetShapes();

    // Loop over shapes
    for (size_t s = 0; s < shapes.size(); s++) {
        // Loop over faces(polygon)
        size_t index_offset = 0;
        for (size_t i = 0; i < shapes[s].mesh.indices.size(); i++){
            unsigned int index = shapes[s].mesh.indices[i].vertex_index;
            out_indices.push_back(index);
        }

        for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
            size_t fv = size_t(shapes[s].mesh.num_face_vertices[f]);
            if (fv != 3) {
                std::cout << "WARNING:TinyObjReader: The mesh " << inputfile << " faces with " << fv <<" vertices!" << std::endl;
            }

            // Loop over vertices in the face.
            for (size_t v = 0; v < fv; v++) {
                // access to vertex
                tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
                tinyobj::real_t vx = attrib.vertices[3 * size_t(idx.vertex_index) + 0];
                tinyobj::real_t vy = attrib.vertices[3 * size_t(idx.vertex_index) + 1];
                tinyobj::real_t vz = attrib.vertices[3 * size_t(idx.vertex_index) + 2];
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
                else {
                    std::cerr << "ERROR: Mesh does not have Normals!" << std::endl;
                }

                // Check if `texcoord_index` is zero or positive. negative = no
                // texcoord data
                if (idx.texcoord_index >= 0) {
                    tinyobj::real_t tx = attrib.texcoords[2 * size_t(idx.texcoord_index) + 0];
                    // Flip the y coordinates to display the textures correctly
                    tinyobj::real_t ty = 1 - attrib.texcoords[2 * size_t(idx.texcoord_index) + 1];
                    out_tex_coord.push_back(tx);
                    out_tex_coord.push_back(ty);
                }
                else {
                    std::cerr << "ERROR: Mesh does not have texture coordinates!" << std::endl;
                }
            }
            index_offset += fv;
        }
    }
}

inline Vec3 compute_tangent_vector(const Vec3& edge1, const Vec3& edge2, const Vec2& deltaUV1, const Vec2& deltaUV2) {
    Vec3 tangent;
    const Scalar denominator = 1.0f / (deltaUV1.x() * deltaUV2.y() - deltaUV1.y() * deltaUV2.x());
    // First vertex tangent
    tangent.x() = denominator * (deltaUV2.y() * edge1.x() - deltaUV1.y() * edge2.x());
    tangent.y() = denominator * (deltaUV2.y() * edge1.y() - deltaUV1.y() * edge2.y());
    tangent.z() = denominator * (deltaUV2.y() * edge1.z() - deltaUV1.y() * edge2.z());
    return tangent.normalized();
}

inline void computeRenderMeshTangentVectors(RenderMesh& mesh) {
    mesh.tangents.resize(mesh.vertices.size());
    const unsigned int nTriangles = static_cast<unsigned int>(mesh.vertices.size()) / 9;
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
    LoadRenderMeshTinyOBJ(filename, vertices, normals, texcoords, indices);
    computeRenderMeshTangentVectors(*this);
    recenter_mesh(*this, compute_COM_position_PARTICLES(vertices));
}

RenderMesh::RenderMesh(const SimulationMesh& simMesh) {
    for (unsigned int i = 0; i < simMesh.indices.size(); i++) {
        unsigned int index = simMesh.indices[i];
        // Vertices
        vertices.push_back(simMesh.vertices[3*index+0]);
        vertices.push_back(simMesh.vertices[3*index+1]);
        vertices.push_back(simMesh.vertices[3*index+2]);

        // All other fields uninitialized
        normals.push_back(0);
        normals.push_back(0);
        normals.push_back(0);

        tangents.push_back(0);
        tangents.push_back(0);
        tangents.push_back(0);

        texcoords.push_back(0);
        texcoords.push_back(0);
    }
    // Compute normals and tangents
    updateFromSimulationMesh(simMesh);
    recenter_mesh(*this, compute_COM_position_PARTICLES(vertices));
}

void LoadVerticesAndIndicesTinyOBJ(std::string inputfile, std::vector<Scalar>& out_vertices, std::vector<unsigned int>& out_indices) {
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

    out_vertices.resize(attrib.vertices.size());
    for (size_t i = 0; i < attrib.vertices.size(); i++) {
        out_vertices[i] = static_cast<Scalar>(attrib.vertices[i]);
    }

    for (size_t s = 0; s < shapes.size(); s++) {
        const tinyobj::shape_t& shape = shapes[s];
        for (unsigned int j = 0; j < shape.mesh.num_face_vertices.size(); j++)
        {
            const int fv = shape.mesh.num_face_vertices[j];
            if (fv != 3) {
                std::cout << "WARNING:TinyObjReader: The mesh " << inputfile << " faces with " << fv <<" vertices!" << std::endl;
            }
        }
        for (size_t i = 0; i < shape.mesh.indices.size(); i++){
            unsigned int index = shape.mesh.indices[i].vertex_index;
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

void recenter_mesh(RenderMesh& mesh, const Vec3& com) {
    for (unsigned int i = 0; i < mesh.vertices.size() / 3; i++) {
        mesh.vertices[3*i+0] -= com.x();
        mesh.vertices[3*i+1] -= com.y();
        mesh.vertices[3*i+2] -= com.z();
    }
}

void RenderMesh::updateFromSimulationMesh(const SimulationMesh& sim_mesh) {
    // The number of triangles of both meshes must be the same
    assert(sim_mesh.indices.size() / 3 == vertices.size() / 9);
    if (tangents.size() == 0)
        tangents.resize(vertices.size());

    indices = sim_mesh.indices;
    const unsigned int nTriangles = static_cast<unsigned int>(sim_mesh.indices.size()) / 3;
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

Vec3 compute_triangle_angles(const Vec3& p1, const Vec3& p2, const Vec3& p3) {const Vec3 AB = p2 - p1;
    const Vec3 AC = p3 - p1;
    const Vec3 BC = p3 - p2;

    const Scalar AC_len = AC.norm();
    const Scalar AB_len = AB.norm();
    const Scalar BC_len = BC.norm();
    const Scalar threshold = 1e-8;
    if (AC_len < threshold or AB_len < threshold or BC_len < threshold) return Vec3::Zero();

    const Scalar cos1 = AB.dot(AC) / (AC_len * AB_len);
    const Scalar cos2 = - BC.dot(AC) / (AC_len * BC_len);

    const Scalar TRIANGLE_TOTAL_ANGLE = M_PI;
    const Scalar angle1 = std::acos(cos1);
    const Scalar angle2 = std::acos(cos2);
    const Scalar angle3 = M_PI - angle1 - angle2;
    return Vec3(angle1, angle2, angle3);
}

// #define ANGLE_WEIGHTED_NORMAL_SMOOTHING
void RenderMesh::smoothNormals() {
    std::vector<Scalar> new_normals;
    new_normals.resize(3*indices.size(), 0.0);

    // Compute smooth normals
    const unsigned int nTriangles = static_cast<unsigned int>(indices.size()) / 3;
    for (unsigned int i = 0; i < nTriangles; i++) {
        // The triangle only has one normal!
        const Vec3 normal = Vec3(normals[9*i + 0 + 3*0], normals[9*i + 1 + 3*0], normals[9*i + 2 + 3*0]);

#ifdef ANGLE_WEIGHTED_NORMAL_SMOOTHING
        const Vec3 vertexA = Vec3(vertices[9*i + 0 + 3*0], vertices[9*i + 1 + 3*0], vertices[9*i + 2 + 3*0]);
        const Vec3 vertexB = Vec3(vertices[9*i + 0 + 3*1], vertices[9*i + 1 + 3*1], vertices[9*i + 2 + 3*1]);
        const Vec3 vertexC = Vec3(vertices[9*i + 0 + 3*2], vertices[9*i + 1 + 3*2], vertices[9*i + 2 + 3*2]);

        const Scalar TRIANGLE_TOTAL_ANGLE = M_PI;
        const Vec3 angles = compute_triangle_angles(vertexA, vertexB, vertexC);
#endif

        for (unsigned int j = 0; j < 3; j++) { // 3 vertices in a triangle
#ifdef ANGLE_WEIGHTED_NORMAL_SMOOTHING
            const Scalar coeff = angles[j] / TRIANGLE_TOTAL_ANGLE;
#else
            const Scalar coeff = 1;
#endif
            const unsigned int idx = 3*indices[3*i+j];
            new_normals[idx + 0] += coeff * normal.x();
            new_normals[idx + 1] += coeff * normal.y();
            new_normals[idx + 2] += coeff * normal.z();
        }
    }

    // Update the normals
    for (unsigned int i = 0; i < nTriangles; i++) {
        for (unsigned int j=0; j < 3; j++) { // j is x, y, z
            // The same normal for each vertex of the triangle
            normals[9*i + j + 3*0] = new_normals[3*indices[3*i + 0] + j];
            normals[9*i + j + 3*1] = new_normals[3*indices[3*i + 1] + j];
            normals[9*i + j + 3*2] = new_normals[3*indices[3*i + 2] + j];
        }
    }
}

std::vector<Scalar> LoadCurveTinyObj(std::string inputfile) {
    std::vector<Scalar> vertices;

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

    vertices.resize(attrib.vertices.size());
    for (size_t s = 0; s < shapes.size(); s++) {
        const tinyobj::shape_t& shape = shapes[s];
        for (size_t i = 0; i < shape.lines.indices.size(); i++){
            unsigned int index = shape.lines.indices[i].vertex_index;
            vertices[3*index] = static_cast<Scalar>(attrib.vertices[3*index]);
            vertices[3*index+1] = static_cast<Scalar>(attrib.vertices[3*index+1]);
            vertices[3*index+2] = static_cast<Scalar>(attrib.vertices[3*index+2]);
        }
    }
    return vertices;
}
