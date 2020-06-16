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
  LATCH = 5,
  MODULE = 6,
};

// This abstractly represents an vertex in the circuit.
class Vertex {
 public:
  Vertex(size_t index, VertexType type)
      : index_(index),
        partition_(-1),
        type_(type) {
      name_ = "not-yet-generated";
  }

  void GenerateName() {
    switch (type_) {
      case VertexType::IN_PIN:
        name_ = "IN_";
        break;
      case VertexType::OUT_PIN:
        name_ = "OUT_";
        break;
      case VertexType::LUT:
        name_ = "LUT_";
        break;
      case VertexType::FLIP_FLOP:
        name_ = "FF_";
        break;
      case VertexType::LATCH:
        name_ = "LATCH_";
        break;
      case VertexType::MODULE:
        name_ = instance_of_ + "_";
        break;
      default:
        name_ = "UNKNOWN_";
        break;
    }
    name_ += std::to_string(index_);
  }

  bool IsSynchronous() const {
    switch (type_) {
      case VertexType::FLIP_FLOP:
      case VertexType::LATCH:
        return true;
      default:
        return false;
    }
  }

  const std::string &name() const { return name_; }
  void set_name(const std::string &name) { name_ = name; }

  size_t index() const { return index_; }
  void set_index(const size_t index) { index_ = index; }

  int partition() const { return partition_; }
  void set_partition(const int partition) { partition_ = partition; }

  double Cost() const { return 1.0; }

  std::vector<Edge*> &in() { return in_; }
  const std::vector<Edge*> &in() const { return in_; }

  std::vector<std::string> &in_ports() { return in_ports_; }
  const std::vector<std::string> &in_ports() const { return in_ports_; }

  std::vector<std::string> &out_ports() { return out_ports_; }
  const std::vector<std::string> &out_ports() const { return out_ports_; }

  std::vector<Edge*> &out() { return out_; }
  const std::vector<Edge*> &out() const { return out_; }

  double weight() { return weight_; }
  void set_weight(const double weight) { weight_ = weight; }

  VertexType type() const { return type_; }
  void set_type(const VertexType type) { type_ = type; }

  const std::string &instance_of() const { return instance_of_; }
  void set_instance_of(const std::string &instance_of) {
    instance_of_ = instance_of;
  }

  const std::string &original_cell_name() const { return original_cell_name_; }
  void set_original_cell_name(const std::string &original_cell_name) {
    original_cell_name_ = original_cell_name;
  }

 private:
  std::string name_;
  size_t index_;
  int partition_;  // -1 means no partition.
  std::vector<Edge*> in_;
  std::vector<Edge*> out_;
  std::vector<std::string> in_ports_;
  std::vector<std::string> out_ports_;
  double weight_;
  VertexType type_;
  std::string instance_of_;
  std::string original_cell_name_;
};

} // namespace jsondotrulo

#endif // _VERTEX_H_
