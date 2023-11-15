#include <cstddef>
#include <vector>
#include <iostream>
#include "raylib.h"
#include "utility_functions.hpp"

#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* .cc
// Optional. define TINYOBJLOADER_USE_MAPBOX_EARCUT gives robust trinagulation. Requires C++11
// #define TINYOBJLOADER_USE_MAPBOX_EARCUT
#include "tiny_obj_loader.h"

Mesh LoadMeshFromVectors(std::vector<unsigned short>& indices, std::vector<float>& vertices, std::vector<float>& normals, std::vector<float>& texcoords) {
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

    std::cout << attrib.vertices.size() << std::endl;
    std::cout << attrib.texcoords.size() << std::endl;

    texcoords.resize(attrib.vertices.size() / 3 * 2, 0);
    normals.resize(attrib.vertices.size());

    for (size_t s = 0; s < shapes.size(); s++) {
        const tinyobj::shape_t& shape = shapes[s];
        for (size_t i = 0; i < shape.mesh.indices.size(); i++){
            size_t index = shape.mesh.indices[i].vertex_index;
            indices.push_back(index);

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
