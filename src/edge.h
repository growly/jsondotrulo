#ifndef _EDGE_H_
#define _EDGE_H_

#include <set>

namespace jsondotrulo {

class Vertex;

// This abstractly represents an edge in the circuit. It's a net.
struct Edge {
  ~Edge() {
    if (name == "fpu_mul_31.421")
      std::cout << "wtf?" << std::endl;
  }
  std::string name;
  double weight;
  size_t index;
  int partition;  // -1 means no partition.
  std::set<Vertex*> in;
  std::set<Vertex*> out;
};

} // namespace jsondotrulo

#endif // _EDGE_H_
