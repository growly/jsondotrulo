#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "edge.h"
#include "path.h"
#include "vertex.h"

namespace jsondotrulo {

// This abstractly represents a self-contained graph, i.e. a model.
// TODO(aryap): This owns its child objects and must delete them.
class Graph {
 public:
  Graph(const std::string &name)
      : name_(name) {}
  ~Graph();

  void AddInputEdges(
      const std::unordered_map<std::string, std::vector<std::string>>
          &input_edges);
  void AddOutputEdges(
      const std::unordered_map<std::string, std::vector<std::string>>
          &output_edges);

  void AddVertex(
        const VertexType &type,
        const std::unordered_map<std::string, std::vector<std::string>>
            &input_edges,
        const std::unordered_map<std::string, std::vector<std::string>>
            &output_edges,
        const std::string &instance_of,
        const std::string &original_cell_name);

  void ExpandInstances(
      const std::unordered_map<std::string, Graph*> &modules_by_name);

  // Finds all paths between synchronous elements (flip flops). Assigns weights
  // accordingly, I guess. TODO(aryap): Find the right weight-assignment
  // scheme.
  void WeightCombinatorialPaths();

  // More generally we might want to do some sanity checking. In this example,
  // we find LUTs with no inputs or with no outputs.
  void FindStrangeLUTs() const;

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

  const std::vector<Vertex*> vertices() const { return vertices_; }

  const std::string &name() const { return name_; }

 private:
  Edge *FindOrCreateEdge(const std::string &name);

  void UpdateEdgeWeightsForPath(
      Path *start_path,
      std::set<Path*> *critical_paths,
      std::unordered_map<Edge*, Path*> *critical_path_by_edge,
      bool *true_on_update);

  std::string name_;
  std::unordered_map<std::string, std::vector<Edge*>> inputs_;
  std::unordered_map<std::string, std::vector<Edge*>> outputs_;

  // The vertex's index in this array is its unique identifier.
  std::vector<Vertex*> vertices_;
  // All of the Edges, in order.
  std::vector<Edge*> edges_;
  // All of the Edges, indexed by their name.
  std::unordered_map<std::string, Edge*> edges_by_name_;
  // All of the Vertices, indexed by their name.
  std::unordered_map<std::string, Vertex*> vertices_by_name_;
  // Edges by partition: if all vertices in/out of an edge are in a single
  // partition, it is listed here. The '-1' partition means 'unpartitioned'.
  std::map<int, std::set<Edge*>> edges_by_partition_;

  // Instances of other modules; these should be replaced with the contents of
  // those modules.
  std::vector<Vertex*> instances_;
};


} // namespace jsondotrulo

#endif // _GRAPH_H_
