#include <assert.h>
#include <iostream>

#include "edge.h"
#include "path.h"
#include "vertex.h"

namespace jsondotrulo {

// We don't enforce loop-free paths here, but the user can by checking
// ContainsVertex and ContainsEdge.
void Path::Append(const std::pair<Vertex*, Edge*> &hop) {
  assert(hop.first != nullptr);
  hops_.push_back(hop);
}

void Path::Append(const Path &path) {
  for (const auto &hop : path) {
    Append(hop);
  }
}

double Path::Cost() const {
  double cost = hops_.size();
  if (source_) cost += source_->Cost();
  return cost;
}

bool Path::PenultimateHop(std::pair<Vertex*, Edge*> *hop) {
  if (hops_.size() < 2) return false;
  *hop = hops_[hops_.size()  - 2];
  return true;
}

bool Path::ContainsVertex(Vertex *vertex) const {
  assert(vertex != nullptr);
  if (source_ && source_->ContainsVertex(vertex)) return true;
  for (const auto &hop : hops_)
    if (hop.first == vertex)
      return true;
  return false;
}

bool Path::ContainsName(const std::string &name) const {
  if (name.empty()) return false;
  if (source_ && source_->ContainsName(name)) return true;
  for (const auto &hop : hops_)
    if (hop.first->name() == name || hop.second->name == name)
      return true;
  return false;
}

bool Path::ContainsEdge(Edge *edge) const {
  assert(edge != nullptr);
  if (source_ && source_->ContainsEdge(edge)) return true;
  for (const auto &hop : hops_)
    if (hop.second == edge)
      return true;
  return false;
}

std::string Path::AsString() const {
  return AsString(true);
}

std::string Path::AsString(bool print_final) const {
  std::string repr;
  if (source_) repr += source_->AsString(false);
  for (auto &pair : hops_) {
    if (pair.second != nullptr)
      repr += pair.first->name() + " -> " + pair.second->name + " -> ";
    else if (print_final)
      repr += pair.first->name();
  }
  return repr;
}

void Path::Print() const {
  std::cout << AsString() << std::endl;
}

} // namespace jsondotrulo
