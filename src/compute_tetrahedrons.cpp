#include "mesh.hpp"
#include <vector>
// #define TETLIBRARY
#include <tetgen.h>


void tetgen_compute_tetrahedrons(const std::vector<unsigned int>& triangle_indices, const std::vector<Scalar>& triangle_vertices, TetrahedronMesh& tmesh) {
    tetgen_compute_tetrahedrons(triangle_indices, triangle_vertices, tmesh.indices, tmesh.vertices);
}

void tetgen_compute_tetrahedrons(const std::vector<unsigned int>& triangle_indices, const std::vector<Scalar>& triangle_vertices,
                                 std::vector<unsigned int>& out_tetrahedron_indices, std::vector<Scalar>& out_tetrahedron_vertices) {

    // Create the tetgen input object
    tetgenio in, out;
    tetgenio::facet *facet;
    tetgenio::polygon *polygon;

    in.firstnumber = 0; // indices start at 0

    // Tetgen deals with REAL == double
    std::vector<double> doubleVertices(triangle_vertices.begin(), triangle_vertices.end());

    // Initialize the triangle mesh vertices
    in.numberofpoints = static_cast<int>(doubleVertices.size()) / 3;
    in.pointlist = new REAL[in.numberofpoints * 3];
    memcpy(in.pointlist, doubleVertices.data(), doubleVertices.size()*sizeof(double));

    // Initialize the triangle mesh indices
    in.numberoffacets = static_cast<int>(triangle_indices.size()) / 3;
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];
    for (int i = 0; i < in.numberoffacets; i++) {
        in.facetmarkerlist[i] = 0; // Initialize this list as zero to treat all the faces equally

        facet = &in.facetlist[i];
        facet->numberofpolygons = 1; // One triangle
        facet->polygonlist = new tetgenio::polygon[facet->numberofpolygons];
        facet->numberofholes = 0;
        facet->holelist = nullptr;

        polygon = &facet->polygonlist[0]; // The triangle in the facet
        polygon->numberofvertices = 3; // The polygon is a triangle
        polygon->vertexlist = new int[polygon->numberofvertices];
        polygon->vertexlist[0] = triangle_indices[3*i+0];
        polygon->vertexlist[1] = triangle_indices[3*i+1];
        polygon->vertexlist[2] = triangle_indices[3*i+2];
    }

    // Generate the tetrahedron mesh
    tetgenbehavior behavior;
    behavior.plc = 1;
    behavior.quality = 1;
    behavior.nobisect = 1; // preserve Surface mesh
    // behavior.minratio = 1.414;
    tetrahedralize(&behavior, &in, &out);

    // Extract the tetrahedra information from "out"
    out_tetrahedron_indices.clear();
    out_tetrahedron_vertices.clear();

    // Tetrahedron indices
    for (int i = 0; i < out.numberoftetrahedra; i++) {
        for (int j = 0; j < 4; j++) {
            out_tetrahedron_indices.push_back(out.tetrahedronlist[i*4+j]);
        }
    }

    // Tetrahedron vertices
    for (int i = 0; i < out.numberofpoints; i++) {
        out_tetrahedron_vertices.push_back(out.pointlist[i*3+0]);
        out_tetrahedron_vertices.push_back(out.pointlist[i*3+1]);
        out_tetrahedron_vertices.push_back(out.pointlist[i*3+2]);
    }
}

TetrahedronMesh::TetrahedronMesh(const SimulationMesh& simMesh) {
    tetgen_compute_tetrahedrons(simMesh.indices,simMesh.vertices, *this);
}
