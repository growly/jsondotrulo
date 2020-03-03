#include "graph.h"
#include <sstream>

#include <iostream>

#include "edge.h"
#include "vertex.h"

namespace jsondotrulo {

Graph::~Graph() {
  for (const Vertex *vertex : vertices_)
    delete vertex;
  for (const Edge *edge : edges_)
    delete edge;
}

Edge *Graph::FindOrCreateEdge(const std::string &name) {
  Edge *net;
  auto edge_it = edges_by_name_.find(name);
  if (edge_it != edges_by_name_.end()) {
    net = edge_it->second;
  } else {
    size_t index = edges_.size();
    net = new Edge({name, index, -1, {}, {}});
    edges_.push_back(net);
    edges_by_name_.insert({name, net});
    edges_by_partition_[-1].insert(net);
  }
  return net;
}

void Graph::AddVertex(
      const std::vector<std::string> &input_edge_names,
      const std::string &output_edge_name) {
  size_t vertex_index = vertices_.size();
  Vertex *vertex = new Vertex(
      {std::to_string(vertex_index), vertex_index, -1, {}, nullptr, 0.0});
  Edge *out_edge = FindOrCreateEdge(output_edge_name);
  out_edge->in.insert(vertex);
  vertex->out = out_edge;
  for (const auto &input : input_edge_names) {
    Edge *in_edge = FindOrCreateEdge(input);
    in_edge->out.insert(vertex);
    vertex->in.push_back(in_edge);
  }
  vertices_.push_back(vertex);
}

void Graph::AddInputEdges(const std::vector<std::string> &input_edge_names) {
  for (const std::string &input : input_edge_names) {
    Edge *net = FindOrCreateEdge(input);
    size_t index = vertices_.size();
    Vertex *in_vertex = new Vertex(
        {"IN_" + std::to_string(index), index, -1, {}, nullptr, 0.0});
    vertices_.push_back(in_vertex);
    net->in.insert(in_vertex);
    inputs_.push_back(net);
  }
}

void Graph::AddOutputEdges(const std::vector<std::string> &output_edge_names) {
  for (const std::string &output : output_edge_names) {
    Edge *net = FindOrCreateEdge(output);
    size_t index = vertices_.size();
    Vertex *out_vertex = new Vertex(
        {"OUT_" + std::to_string(index), index, -1, {}, nullptr, 0.0});
    vertices_.push_back(out_vertex);
    net->out.insert(out_vertex);
    outputs_.push_back(net);
  }
}

void Graph::ReadHMETISPartitions(const std::string &in_str) {
  std::istringstream input(in_str);
  size_t vertex_index = 0;
  for (std::string line; std::getline(input, line); ) {
    int partition = std::stoi(line);
    vertices_[vertex_index]->partition = partition;
    ++vertex_index;
  }

  // Now that all partitions are known, assign edges to the right partitions.
  edges_by_partition_.clear();
  for (Edge *edge : edges_) {
    int partition = -1;
    if (!edge->in.empty()) {
      partition = (*edge->in.begin())->partition;
    } else if (!edge->out.empty()) {
      partition = (*edge->out.begin())->partition;
    }

    for (const Vertex *vertex : edge->in) {
      if (vertex->partition != partition) {
        partition = -1;
      }
    }
    for (const Vertex *vertex : edge->out) {
      if (vertex->partition != partition) {
        partition = -1;
      }
    }
    edge->partition = partition;
    edges_by_partition_[partition].insert(edge);
  }
}

void Graph::Print() const {
  std::cout << "Graph \"" << name_ << "\":" << std::endl
      << "\tinputs: " << inputs_.size() << std::endl;
  for (const auto &edge : inputs_)
    std::cout << "\t\t" << edge->name << std::endl;
  std::cout << "\toutputs: " << outputs_.size() << std::endl;
  for (const auto &edge : outputs_)
    std::cout << "\t\t" << edge->name << std::endl;

  std::cout << "\tvertices: " << vertices_.size() << std::endl;
  std::cout << "\tedges: " << edges_.size() << std::endl;
  std::cout << "\tpartitions: "
            << std::to_string(edges_by_partition_.size() - 1) << std::endl;
}

std::string Graph::AsDOT() const {
  std::string out_str = "digraph " + name_ + " {\n";
  for (const auto &item : edges_by_partition_) {
    int partition = item.first;
    if (partition >= 0) {
      // You actually have to use the keyword "cluster" as a prefix:
      out_str += "subgraph cluster_" + std::to_string(partition) +
                 " {\n  label=\"partition_" + std::to_string(partition) +
                 "\";\n  node [style=filled];\n  color=\"black\";\n";
    }
    for (const Edge *edge : item.second) {
      for (const Vertex *in : edge->in) {
        for (const Vertex *out : edge->out) {
          out_str += "  " + in->name + " -> " + out->name + " [label=\""
                     + edge->name + "\"];\n";
        }
      }
    }
    if (partition >= 0) {
      out_str += "}\n";
    }
  }
  out_str += "}\n";
  return out_str;
}

std::string Graph::AsMFile() const {
  size_t num_vertices = vertices_.size();
  //int **matrix = new int*[num_vertices];
  //for (size_t i = 0; i < num_vertices; ++i) {
  //  matrix[i] = new int[num_vertices];
  //}
  // Some C++11 magic?
  int matrix[num_vertices][num_vertices] = {0};

  for (const Edge *edge : edges_) {
    for (const Vertex *in : edge->in) {
      for (const Vertex *out : edge->out) {
        matrix[in->index][out->index] = 1;
        matrix[out->index][in->index] = 1;
      }
    }
  }

  std::string out_str = name_ + "_adjacency_matrix = [";
  for (size_t i = 0; i < num_vertices; ++i) {
    for (size_t j = 0; j < num_vertices; ++j) {
      out_str += std::to_string(matrix[i][j]);
      if (j < num_vertices - 1) out_str += " ";
    }
    out_str += ";\n";
  }
  out_str += "]\n";
  return out_str;
}

std::string Graph::AsHMETIS() const {
  std::string out_str = "\% num_hyperedges num_vertices fmt\n";
  out_str += std::to_string(edges_.size()) + " "
             + std::to_string(vertices_.size()) + "\n";

  // From hMETIS manual:
  // The remaining |E| lines store the vertices contained in each hyperedge -
  // one line per hyperedge. In particular, the ith line (excluding comment
  // lines) contains the vertices that are included in the (i-1)th hyperedge.
  //
  // These hyperedges are treated as UNDIRECTED.
  out_str += "\% vertices connected by hyperedges, one line per hyperedge\n";
  for (const Edge *edge : edges_) {
    for (Vertex *vertex : edge->in)
      out_str += std::to_string(vertex->index + 1) + " ";
    for (Vertex *vertex : edge->out)
      out_str += std::to_string(vertex->index + 1) + " ";
    out_str[out_str.size() - 1] = '\n';
  }
  return out_str;
}

} // namespace jsondotrulo
