#ifndef EDGE_H_
#define EDGE_H_

#include <unordered_map>

class Edge{
    public:
        int a, b;

        int opposite;

        Edge() :a(-1) ,b(-1), opposite(-1) {};
        Edge(int a, int b, int opposite)
            : a(a), b(b), opposite(opposite) {}

        Edge reversed() { return Edge(b, a, -1); }
};

// Have a way to hash the Edge class
namespace std{
  template<>
  struct hash<Edge>{
    unsigned int operator()(const Edge& key) const;
  };
}

bool operator==(const Edge& e1, const Edge& e2);

#endif // EDGE_H_
