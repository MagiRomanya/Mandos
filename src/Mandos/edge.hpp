#ifndef MANDOS_EDGE_H_
#define MANDOS_EDGE_H_

#include <unordered_map>

namespace mandos
{
class Edge
{
public:
    int a, b;

    int opposite;

    Edge()
        : a(-1)
        , b(-1)
        , opposite(-1){};
    Edge(int a, int b, int opposite)
        : a(a)
        , b(b)
        , opposite(opposite)
    {
    }

    Edge reversed()
    {
        return Edge(b, a, -1);
    }
};
}  // namespace mandos

// Have a way to hash the Edge class
// TODO This might be UB
namespace std{
template <>
struct hash<mandos::Edge> {
    unsigned int operator()(const mandos::Edge& key) const;
};
}

bool operator==(const mandos::Edge& e1, const mandos::Edge& e2);

#endif  // MANDOS_EDGE_H_
