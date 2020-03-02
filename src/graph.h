#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <string>
#include <unordered_map>
#include <vector>

namespace jsondotrulo {

class Edge;
class Vertex;

// This abstractly represents a self-contained graph, i.e. a model.
// TODO(aryap): This owns its child objects and must delete them.
class Graph {
 public:
  Graph(const std::string &name)
      : name_(name) {}
  ~Graph();

  void AddInputEdges(const std::vector<std::string> &input_edge_names);
  void AddOutputEdges(const std::vector<std::string> &output_edge_names);

  void AddVertex(
      const std::vector<std::string> &input_edge_names,
      const std::string &output_edge_name);

  void Print() const;

  std::string AsDOT() const;

  std::string AsMFile() const;

  const std::string &name() const { return name_; }

 private:
  Edge *FindOrCreateEdge(const std::string &name);

  std::string name_;
  std::vector<Edge*> inputs_;
  std::vector<Edge*> outputs_;

  // The vertex's index in this array is its unique identifier.
  std::vector<Vertex*> vertices_;
  std::unordered_map<std::string, Edge*> edges_;
};


} // namespace jsondotrulo

#endif // _GRAPH_H_
