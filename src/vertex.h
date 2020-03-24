#ifndef _VERTEX_H_
#define _VERTEX_H_

#include <string>

namespace jsondotrulo {

struct Edge;

enum VertexType {
  UNKNOWN = 0,
  IN_PIN = 1,
  OUT_PIN = 2,
  LUT = 3,
  FLIP_FLOP = 4,
};

// This abstractly represents an vertex in the circuit.
class Vertex {
 public:
  Vertex(size_t index, VertexType type)
      : index_(index),
        partition_(-1),
        type_(type) {
      name_ = GenerateName();
  }

  const std::string &name() const { return name_; }
  void set_name(const std::string &name) { name_ = name; }

  size_t index() const { return index_; }
  void set_index(const size_t index) { index_ = index; }

  int partition() const { return partition_; }
  void set_partition(const int partition) { partition_ = partition; }

  std::vector<Edge*> &in() { return in_; }
  const std::vector<Edge*> &in() const { return in_; }

  Edge *out() { return out_; }
  const Edge *out() const { return out_; }
  void set_out(Edge *out) { out_ = out; }

  double weight() { return weight_; }
  void set_weight(const double weight) { weight_ = weight; }

  VertexType type() const { return type_; }
  void set_type(const VertexType type) { type_ = type; }

 private:
  std::string GenerateName() const {
    std::string name;
    switch (type_) {
      case VertexType::IN_PIN:
        name = "IN_";
        break;
      case VertexType::OUT_PIN:
        name = "OUT_";
        break;
      case VertexType::LUT:
        name = "LUT_";
        break;
      case VertexType::FLIP_FLOP:
        name = "FF_";
        break;
      default:
        name = "UNKNOWN_";
        break;
    }
    name += std::to_string(index_);
    return name;
  }

  std::string name_;
  size_t index_;
  int partition_;  // -1 means no partition.
  std::vector<Edge*> in_;
  Edge *out_;
  double weight_;
  VertexType type_;
};

} // namespace jsondotrulo

#endif // _VERTEX_H_
