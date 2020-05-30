#include <iostream>

#include "edge.h"
#include "path.h"
#include "vertex.h"

namespace jsondotrulo {

// We don't enforce loop-free paths here, but the user can by checking
// ContainsVertex and ContainsEdge.
void Path::Append(const std::pair<Vertex*, Edge*> &hop) {
  hops_.push_back(hop);
}

double Path::Cost() const {
  return hops_.size();
}

bool Path::ContainsVertex(Vertex *vertex) const {
  return vertices_.find(vertex) != vertices_.end();
}

bool Path::ContainsEdge(Edge *edge) const {
  return edges_.find(edge) != edges_.end();
}

std::string Path::AsString() const {
  std::string repr;
  for (auto &pair : hops_) {
    repr += pair.first->name();
    if (pair.second != nullptr)
      repr += " -> " + pair.second->name + " -> ";
    else
      repr += " -> null";
  }
  repr += "]";
  return repr;
}

void Path::Print() const {
  std::cout << AsString() << std::endl;
}

} // namespace jsondotrulo