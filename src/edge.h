#ifndef _EDGE_H_
#define _EDGE_H_

#include <set>

namespace jsondotrulo {

class Vertex;

// This abstractly represents an edge in the circuit. It's a net.
struct Edge {
  std::string name;
  std::set<Vertex*> in;
  std::set<Vertex*> out;
};

} // namespace jsondotrulo

#endif // _EDGE_H_
