#include "mesh.hpp"
#include "edge.hpp"
#include "linear_algebra.hpp"
#include "utility_functions.hpp"
#include <cassert>
#include <cstdio>
#include <unordered_set>
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
    std::unordered_map<Edge, int> in_map;

    unsigned int externalSize = 0;
    unsigned int internalSize = 0;

    for (size_t i=0; i < indices.size(); i+=3) {
        const unsigned int a = indices[i];
        const unsigned int b = indices[i+1];
        const unsigned int c = indices[i+2];
        // Edge a, b
        // If the inverse edge has already been added, it adds the current edge next to it
        Edge e = Edge(a, b, c);
        if (in_map.count(e.reversed())){
            int inverse_edge_index = in_map[e.reversed()];
            edges[inverse_edge_index + 1] = e;
            internalSize+=2;
        }
        else{
            in_map[e] = edges.size();
            edges.push_back(e);
            edges.push_back(Edge(-1,-1,-1)); // Add a dummy edge to the list for later removal if it's not overwritten
        }

        // Edge c, a
        e = Edge(c, a, b);
        if (in_map.count(e.reversed())){
            int inverse_edge_index = in_map[e.reversed()];
            edges[inverse_edge_index + 1] = e;
            internalSize+=2;
        }
        else{
            in_map[e] = edges.size();
            edges.push_back(e);
            edges.push_back(Edge(-1,-1,-1));
        }

        // Edge b, c
        e = Edge(b, c, a);
        if (in_map.count(e.reversed())){
            int inverse_edge_index = in_map[e.reversed()];
            edges[inverse_edge_index + 1] = e;
            internalSize+=2;
        }
        else{
            in_map[e] = edges.size();
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


struct Triangle {
    unsigned int a, b, c;
};

bool operator==(const Triangle& t1, const Triangle& t2) {
    return ((t1.a == t2.a) && (t1.b == t2.b) && (t1.c == t2.c)) ||
        ((t1.a == t2.b) && (t1.b == t2.c) && (t1.c == t2.a)) ||
        ((t1.a == t2.c) && (t1.b == t2.a) && (t1.c == t2.b));
}

namespace std {
    template<>
    struct hash<Triangle>{
        unsigned int operator()(const Triangle& key) const {
            size_t hashValue = 17; // Choose a prime number as a seed
            hashValue = hashValue * 31 + key.a;
            hashValue = hashValue * 31 + key.b;
            hashValue = hashValue * 31 + key.c;
            return hashValue;
        }
    };
}

void compute_triangle_indices_from_tetrahedron_indices(const std::vector<unsigned int>& tet_ind, std::vector<unsigned int>& out_ind) {
    assert(tet_ind.size() % 4 == 0);

    // The tetrahderon mesh will have repeated triangles
    std::unordered_set<Triangle> known_triangles;

    // Iterate tetrahedrons
    for (unsigned int i = 0; i < tet_ind.size()/4; i++) {
        const unsigned int a = tet_ind[4*i+1];
        const unsigned int b = tet_ind[4*i+0];
        const unsigned int c = tet_ind[4*i+2];
        const unsigned int d = tet_ind[4*i+3];

        // Each tetrahedron has 4 triangular faces
        const Triangle t1 = Triangle(a,b,c);
        const Triangle t2 = Triangle(a,d,b);
        const Triangle t3 = Triangle(a,c,d);
        const Triangle t4 = Triangle(b,d,c);

        // Add the triangle to the output when it is not repeated
        const Triangle triangles[] = {t1,t2,t3,t4};
        for (unsigned int j=0; j < 4; j++) {
            const Triangle& t = triangles[j];
            if (true or known_triangles.contains(t)) {
                known_triangles.insert(t);
                out_ind.push_back(t.a);
                out_ind.push_back(t.b);
                out_ind.push_back(t.c);
            }
        }
    }
}
