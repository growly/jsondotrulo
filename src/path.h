#ifndef _PATH_H_
#define _PATH_H_

#include <memory>
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

  Path(const std::pair<Vertex*, Edge*> start)
    : source_(nullptr),
      num_edges_for_which_critical_(0) {
    Append(start);
  }

  Path(std::shared_ptr<Path> source)
    : source_(source),
      num_edges_for_which_critical_(0) {}

  Path(const Path &other)
      : hops_(other.hops_.begin(), other.hops_.end()),
        source_(other.source_),
        num_edges_for_which_critical_(0) {
  }

  ~Path() {
    //std::cout << "path deleted: " << AsString() << std::endl;
  }

  void Append(const std::pair<Vertex*, Edge*> &hop);

  void Append(const Path &path);

  double Cost() const;

  std::string AsString() const;
  std::string AsString(bool print_final) const;

  void Print() const;

  size_t size() const {
    return hops_.size();
  }

  bool empty() const {
    return hops_.empty();
  }

  std::pair<Vertex*, Edge*> &front() { return hops_.front(); }
  const std::pair<Vertex*, Edge*> &front() const { return hops_.front(); }

  std::vector<std::pair<Vertex*, Edge*>>::iterator begin() {
    return hops_.begin();
  }
  std::vector<std::pair<Vertex*, Edge*>>::iterator end() {
    return hops_.end();
  }

  std::vector<std::pair<Vertex*, Edge*>>::const_iterator begin() const {
    return hops_.begin();
  }
  std::vector<std::pair<Vertex*, Edge*>>::const_iterator end() const {
    return hops_.end();
  }

  std::shared_ptr<Path> &source() { return source_; }

  // Places a copy of the second-to-last hop in *hop and return true, unless
  // the Path isn't long enough to have a second-to-last hop, in which case
  // *hop is not modified and this returns false.
  bool PenultimateHop(std::pair<Vertex*, Edge*> *hop);

  bool ContainsVertex(Vertex *vertex) const;
  bool ContainsEdge(Edge *edge) const;

  void IncrementNumEdgesForWhichCritical() {
    num_edges_for_which_critical_++;
  }

  void DecrementNumEdgesForWhichCritical() {
    num_edges_for_which_critical_--;
  }

  size_t num_edges_for_which_critical() {
    return num_edges_for_which_critical_;
  }

 private:
  std::vector<std::pair<Vertex*, Edge*>> hops_;

  // Paths are pretty short, so we can afford an O(n) search to avoid the
  // memory-footprint of a log(n) search structure (such as additional sets).
 
  // Sourth Path, if any. This path is treated as part of this path, so is
  // recursively searched for elements and counts towards this Path's length, if
  // any.
  std::shared_ptr<Path> source_;

  size_t num_edges_for_which_critical_;
};

} // namespace jsondotrulo

#endif // _PATH_H_
