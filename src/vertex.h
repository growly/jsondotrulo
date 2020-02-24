#ifndef _VERTEX_H_
#define _VERTEX_H_

namespace jsondotrulo {

struct Edge;

// This abstractly represents an vertex in the circuit.
struct Vertex {
  std::string name;
  std::vector<Edge*> in;
  Edge* out;
  double weight;
};

} // namespace jsondotrulo

#endif // _VERTEX_H_
