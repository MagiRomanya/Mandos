#include <cstddef>
#include <vector>
#include <iostream>
#include "raylib.h"
#include "utility_functions.hpp"

#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* .cc
// Optional. define TINYOBJLOADER_USE_MAPBOX_EARCUT gives robust trinagulation. Requires C++11
// #define TINYOBJLOADER_USE_MAPBOX_EARCUT
#include "tiny_obj_loader.h"

Mesh LoadMeshFromVectors(std::vector<unsigned short>& indices, std::vector<float>& vertices, std::vector<float>& normals, std::vector<float>& texcoord) {
    Mesh mesh = {0};
    mesh.vertexCount = vertices.size() / 3;
    mesh.triangleCount = indices.size() / 3;
    mesh.vertices = (float *)std::memcpy(RL_MALLOC(vertices.size() * sizeof(float)), vertices.data(), vertices.size() * sizeof(float));
    mesh.texcoords = (float *)std::memcpy(RL_MALLOC(texcoord.size() * sizeof(float)), texcoord.data(), texcoord.size() * sizeof(float));
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

    // Loop over faces(polygon)
    size_t index_offset = 0;
    size_t s = 0;
    for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
        size_t fv = size_t(shapes[s].mesh.num_face_vertices[f]);

        // Loop over vertices in the face.
        for (size_t v = 0; v < fv; v++) {
            // access to vertex
            tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];

            tinyobj::real_t vx = attrib.vertices[3*size_t(idx.vertex_index)+0];
            tinyobj::real_t vy = attrib.vertices[3*size_t(idx.vertex_index)+1];
            tinyobj::real_t vz = attrib.vertices[3*size_t(idx.vertex_index)+2];

            // Check if `normal_index` is zero or positive. negative = no normal data
            if (idx.normal_index >= 0) {
                tinyobj::real_t nx = attrib.normals[3*size_t(idx.normal_index)+0];
                tinyobj::real_t ny = attrib.normals[3*size_t(idx.normal_index)+1];
                tinyobj::real_t nz = attrib.normals[3*size_t(idx.normal_index)+2];
            }

            // Check if `texcoord_index` is zero or positive. negative = no texcoord data
            if (idx.texcoord_index >= 0) {
                tinyobj::real_t tx = attrib.texcoords[2*size_t(idx.texcoord_index)+0];
                tinyobj::real_t ty = attrib.texcoords[2*size_t(idx.texcoord_index)+1];
            }
        }
        index_offset += fv;
    }
    vertices = attrib.vertices;
    return LoadMeshFromVectors(indices, vertices, normals, texcoords);
}
