#ifndef _PATH_H_
#define _PATH_H_

#include <set>
#include <string>
#include <utility>
#include <vector>

namespace jsondotrulo {

class Vertex;
struct Edge;

class Path {
 public:
  Path() {}
  Path(const std::pair<Vertex*, Edge*> start) {
    Append(start);
  }
  Path(const Path &other)
      : hops_(other.hops_.begin(), other.hops_.end()),
        vertices_(other.vertices_.begin(), other.vertices_.end()),
        edges_(other.edges_.begin(), other.edges_.end()) {
  }
  ~Path() = default;

  void Append(const std::pair<Vertex*, Edge*> &hop);

  double Cost() const;

  std::string AsString() const;

  void Print() const;

  std::vector<std::pair<Vertex*, Edge*>>::iterator begin() {
    return hops_.begin();
  }
  std::vector<std::pair<Vertex*, Edge*>>::iterator end() {
    return hops_.end();
  }

  bool ContainsVertex(Vertex *vertex) const;
  bool ContainsEdge(Edge *edge) const;

 private:
  std::vector<std::pair<Vertex*, Edge*>> hops_;

  // This structures help us prevent loops faster, by reducing the search time
  // when checking for loops.
  std::set<Vertex*> vertices_;
  std::set<Edge*> edges_;
};

} // namespace jsondotrulo

#endif // _PATH_H_
