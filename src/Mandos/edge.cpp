#include <Mandos/edge.hpp>

unsigned int std::hash<mandos::Edge>::operator()(const mandos::Edge& key) const
{
    return 100000 * key.a + key.b;
}

namespace mandos
{

bool operator==(const Edge& e1, const Edge& e2)
{
    return (e1.a == e2.a) and (e1.b == e2.b);
}

}  // namespace mandos