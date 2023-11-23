#include <vector>
#include <tetgen.h>

#include "utility_functions.hpp"

void tetgen_compute_tetrahedrons(const std::vector<unsigned int>& triangle_indices, const std::vector<float>& triangle_vertices,
                                 std::vector<unsigned int>& out_tetrahedron_indices, std::vector<float>& out_tetrahedron_vertices) {

    // Create the tetgen input object
    tetgenio in, out;
    tetgenio::facet *facet;
    tetgenio::polygon *polygon;

    in.firstnumber = 0; // indices start at 0

    // Tetgen deals with REAL == double
    std::vector<double> doubleVertices(triangle_vertices.begin(), triangle_vertices.end());

    // Initialize the triangle mesh vertices
    in.numberofpoints = doubleVertices.size() / 3;
    in.pointlist = new REAL[in.numberofpoints * 3];
    memcpy(in.pointlist, doubleVertices.data(), doubleVertices.size()*sizeof(double));

    // Initialize the triangle mesh indices
    in.numberoffacets = triangle_indices.size() / 3;
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];
    for (int i = 0; i < in.numberoffacets; i++) {
        in.facetmarkerlist[i] = 0; // Initialize this list as zero to treat all the faces equally

        facet = &in.facetlist[i];
        facet->numberofpolygons = 1; // One triangle
        facet->polygonlist = new tetgenio::polygon[facet->numberofpolygons];
        facet->numberofholes = 0;
        facet->holelist = nullptr;

        polygon = &facet->polygonlist[0];
        polygon->numberofvertices = 3; // The polygon is a triangle
        polygon->vertexlist = new int[polygon->numberofvertices];
        polygon->vertexlist[0] = triangle_indices[3*i+0];
        polygon->vertexlist[1] = triangle_indices[3*i+1];
        polygon->vertexlist[2] = triangle_indices[3*i+2];
    }

    // Generate the tetrahedron mesh
    tetgenbehavior behaviour;
    behaviour.plc = 1;
    behaviour.quality = 1;
    // behaviour.minratio = 1.414;
    tetrahedralize(&behaviour, &in, &out);

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
