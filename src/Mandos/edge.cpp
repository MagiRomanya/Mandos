#include "edge.hpp"
// Have a way to hash the Edge class
unsigned int std::hash<Edge>::operator()(const Edge& key) const{
    return 100000 * key.a + key.b;
}

bool operator==(const Edge& e1, const Edge& e2){
    return (e1.a == e2.a) and (e1.b == e2.b);
}
