#include "mesh_boundary.hpp"
#include "edge.hpp"
#include "linear_algebra.hpp"
#include <cstdio>
#include <vector>

void mesh_boundary(const std::vector<Scalar>& vertices,
                   const std::vector<unsigned int>& indices,
                   std::vector<Edge> &internalEdges,
                   std::vector<Edge> &externalEdges) {
    /* Fills the internal edges and external edges vectors.
     * In the internal edges vectors, we order the vector
     * by puting a semi-edge and its inverse next to each other. */

    std::vector<Edge> edges;
    edges.reserve(2 * vertices.size()); // In the worst case scenario all the mesh will be boudary
    std::unordered_map<Edge, int> in_vector;

    unsigned int externalSize = 0;
    unsigned int internalSize = 0;

    for (size_t i=0; i < indices.size(); i+=3) {
        const unsigned int a = indices[i];
        const unsigned int b = indices[i+1];
        const unsigned int c = indices[i+2];
        // Edge a, b
        // If the inverse edge has already been added, it adds the current edge next to it
        Edge e = Edge(a, b, c);
        if (in_vector.count(e.reversed())){
            int inverse_edge_index = in_vector[e.reversed()];
            edges[inverse_edge_index + 1] = e;
            internalSize+=2;
        }
        else{
            in_vector[e] = edges.size();
            edges.push_back(e);
            edges.push_back(Edge(-1,-1,-1)); // Add a dummy edge to the list for later removal if it's not overwritten
        }

        // Edge c, a
        e = Edge(c, a, b);
        if (in_vector.count(e.reversed())){
            int inverse_edge_index = in_vector[e.reversed()];
            edges[inverse_edge_index + 1] = e;
            internalSize+=2;
        }
        else{
            in_vector[e] = edges.size();
            edges.push_back(e);
            edges.push_back(Edge(-1,-1,-1));
        }

        // Edge b, c
        e = Edge(b, c, a);
        if (in_vector.count(e.reversed())){
            int inverse_edge_index = in_vector[e.reversed()];
            edges[inverse_edge_index + 1] = e;
            internalSize+=2;
        }
        else{
            in_vector[e] = edges.size();
            edges.push_back(e);
            edges.push_back(Edge(-1,-1,-1)); // dummy edge
        }
    }

    externalSize = indices.size() - internalSize;

    externalEdges.reserve(externalSize);
    internalEdges.reserve(internalSize);

    // Now we have to remove the dummy edges created by the edge of the mesh where no inverse exist ( the boundary )
    for (size_t i=0; i < edges.size(); i+=2){
        if (edges[i+1].a < 0){ // If the 2nd edge is dummy we have found a external edge with no inverse
            externalEdges.push_back(edges[i]);
        }
        else {
            internalEdges.push_back(edges[i]);
            internalEdges.push_back(edges[i+1]);
        }
    }
}

std::array<unsigned int, 2> count_springs(const std::vector<Scalar>& vertices, const std::vector<unsigned int>& indices) {
    std::vector<Edge> internalEdges, externalEdges;
    mesh_boundary(vertices, indices, internalEdges, externalEdges);

    unsigned int n_flex = internalEdges.size() / 2.0 + externalEdges.size();
    unsigned int n_bend = internalEdges.size() / 2.0;

    return {n_flex, n_bend};
}
