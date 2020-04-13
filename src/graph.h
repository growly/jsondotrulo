#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "edge.h"
#include "vertex.h"

namespace jsondotrulo {

typedef std::vector<std::pair<Vertex*, Edge*>> Path;

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
        const VertexType &type,
        const std::vector<std::string> &input_edge_names,
        const std::vector<std::string> &output_edge_names);

  // Finds all paths between synchronous elements (flip flops). Assigns weights
  // accordingly, I guess. TODO(aryap): Find the right weight-assignment
  // scheme.
  void WeightCombinatorialPaths() ;

  void Print() const;

  std::string AsDOT() const;

  std::string AsMFile() const;

  std::string AsHMETIS() const;

  // TODO(aryap): Verify this actually works.
  std::string AsGraph6() const;

  std::string AsEdgeListWithWeights() const;

  // Read parition information from an hMETIS-format partition file. It is
  // assumed that the hyperedge indices in the file correspond to those read
  // into this graph.
  void ReadHMETISPartitions(const std::string &in_str);

  const std::string &name() const { return name_; }

 private:
  static double ApproximateCost(const Path &path);

  Edge *FindOrCreateEdge(const std::string &name);

  std::string name_;
  std::vector<Edge*> inputs_;
  std::vector<Edge*> outputs_;

  // The vertex's index in this array is its unique identifier.
  std::vector<Vertex*> vertices_;
  // All of the Edges, in order.
  std::vector<Edge*> edges_;
  // All of the Edges, indexed by their name.
  std::unordered_map<std::string, Edge*> edges_by_name_;

  // Edges by partition: if all vertices in/out of an edge are in a single
  // partition, it is listed here. The '-1' partition means 'unpartitioned'.
  std::map<int, std::set<Edge*>> edges_by_partition_;
};


} // namespace jsondotrulo

#endif // _GRAPH_H_
