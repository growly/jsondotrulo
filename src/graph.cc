#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <set>
#include <sstream>
#include <unordered_map>
#include <utility>

#include "graph.h"
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
    net = new Edge({name, 1.0, index, -1, {}, {}});
    edges_.push_back(net);
    edges_by_name_.insert({name, net});
    edges_by_partition_[-1].insert(net);
  }
  return net;
}

double Graph::ApproximateCost(const std::vector<Vertex*> &path) {
  return static_cast<double>(path.size());
}

void Graph::AddVertex(
      const VertexType &type,
      const std::vector<std::string> &input_edge_names,
      const std::string &output_edge_name) {
  size_t vertex_index = vertices_.size();
  Vertex *vertex = new Vertex(vertex_index, type);
  Edge *out_edge = FindOrCreateEdge(output_edge_name);
  out_edge->in.insert(vertex);
  vertex->set_out(out_edge);
  for (const auto &input : input_edge_names) {
    Edge *in_edge = FindOrCreateEdge(input);
    in_edge->out.insert(vertex);
    vertex->in().push_back(in_edge);
  }
  vertices_.push_back(vertex);
}

void Graph::WeightCombinatorialPaths() {
  // The goal is to find every path from one flip to another flip flop, through
  // whatever other vertices.

  // TODO(aryap): Paths are a list of vertices, since for now each vertex has
  // only one output. 'paths' contains every path created in this function. The
  // contents must be deleted at the end. unique_ptr?
  std::vector<std::vector<Vertex*>*> paths;
  std::vector<std::vector<Vertex*>*> complete;
  for (Vertex *start : vertices_) {
    if (start->type() != VertexType::FLIP_FLOP) continue;

    std::set<Vertex*> visited;
    std::vector<Vertex*> *start_path = new std::vector<Vertex*>{start};
    paths.push_back(start_path);
    std::vector<std::pair<Vertex*, std::vector<Vertex*>*>> to_visit;
    to_visit.push_back(std::make_pair(start, start_path));

    // We do not care for the other drivers of the same net. Presumably they
    // are not on the path.
    while (!to_visit.empty()) {
      Vertex *current = to_visit.back().first;
      std::vector<Vertex*> *path = to_visit.back().second;
      to_visit.pop_back();
      if (visited.find(current) != visited.end()) continue;
      visited.insert(current);

      // Paths end at flip flops, but we have to make sure we go somewhere
      // first (otherwise the starting node, which is also a flip flop, would
      // terminate the path). At the same time, loops are ok.
      if (path->size() > 1 && current->type() == VertexType::FLIP_FLOP) {
        // Paths end here.
        complete.push_back(path);
        continue;
      }

      if (current->out() == nullptr) continue;

      for (Vertex *next_hop : current->out()->out) {
        // Create a new path for the each next-hop.
        std::vector<Vertex*> *new_path =
            new std::vector<Vertex*>(path->begin(), path->end());
        new_path->push_back(next_hop);
        paths.push_back(new_path);
        to_visit.push_back(std::make_pair(next_hop, new_path));
      }
    }
  }

  for (std::vector<Vertex*> *path : complete) {
    //std::cout << "cb path: ";
    //for (int i = 0; i < path->size() - 1; ++i) {
    //  std::cout << (*path)[i]->name() << " -> ";
    //}
    //std::cout << path->back()->name() << std::endl;

    // Last hop in path is terminal
    double path_cost = ApproximateCost(*path);
    for (int i = 0; i < path->size() - 1; ++i) {
      Edge *edge = path->at(i)->out();
      if (edge == nullptr) continue;
      edge->weight = std::max(edge->weight, path_cost);
    }
    delete path;
  }
}

void Graph::AddInputEdges(const std::vector<std::string> &input_edge_names) {
  for (const std::string &input : input_edge_names) {
    Edge *net = FindOrCreateEdge(input);
    size_t index = vertices_.size();
    Vertex *in_vertex = new Vertex(index, VertexType::IN_PIN);
    vertices_.push_back(in_vertex);
    net->in.insert(in_vertex);
    inputs_.push_back(net);
  }
}

void Graph::AddOutputEdges(const std::vector<std::string> &output_edge_names) {
  for (const std::string &output : output_edge_names) {
    Edge *net = FindOrCreateEdge(output);
    size_t index = vertices_.size();
    Vertex *out_vertex = new Vertex(index, VertexType::OUT_PIN);
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
    vertices_[vertex_index]->set_partition(partition);
    ++vertex_index;
  }

  // Now that all partitions are known, assign edges to the right partitions.
  edges_by_partition_.clear();
  for (Edge *edge : edges_) {
    int partition = -1;
    if (!edge->in.empty()) {
      partition = (*edge->in.begin())->partition();
    } else if (!edge->out.empty()) {
      partition = (*edge->out.begin())->partition();
    }

    for (const Vertex *vertex : edge->in) {
      if (vertex->partition() != partition) {
        partition = -1;
      }
    }
    for (const Vertex *vertex : edge->out) {
      if (vertex->partition() != partition) {
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
  // dot graph names can't have '/' in them:
  std::replace(out_str.begin(), out_str.end(), '/', '_');
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
          out_str += "  " + in->name() + " -> " + out->name() + " [label=\""
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
        matrix[in->index()][out->index()] = 1;
        matrix[out->index()][in->index()] = 1;
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
      out_str += std::to_string(vertex->index() + 1) + " ";
    for (Vertex *vertex : edge->out)
      out_str += std::to_string(vertex->index() + 1) + " ";
    out_str[out_str.size() - 1] = '\n';
  }
  return out_str;
}

namespace {

// Implements the R(x) function of a bitvector, x, defined for the Graph6
// format here:
//    http://users.cecs.anu.edu.au/~bdm/data/formats.txt
static void Graph6Rx(std::vector<bool> bits, std::vector<uint8_t> *representation) {
  double k = bits.size();
  std::cout << "input bitvector length k = " << k << std::endl;
  size_t num_bytes = static_cast<size_t>(std::ceil(k/6.0));
  std::cout << "num_bytes = " << num_bytes << std::endl;

  // Pad out the input vector with 0s at the right.
  for (int i = k; i < num_bytes; ++i)
    bits.push_back(false);

  // Split into groups of 6 bits each.
  for (int i = 0; i < num_bytes; ++i) {
    uint8_t byte = 0;
    for (int j = 0; j < 6; ++j) {
      int m = i * 6 + j;
      if (bits[m]) {
        byte |= (1U << (5 - j));
      }
    }
    // Add 63 to the byte for this group.
    byte += 63;
    representation->push_back(byte);
  }

  // Print.
  //for (int i = 0; i < representation->size(); ++i) {
  //  std::cout << std::to_string(representation->at(i)) << " ";
  //}
  //std::cout << std::endl;
}

// TODO(aryap): LOL TEST THIS
static void IntegerToBits(uint64_t integer, std::vector<bool> *bits) {
  // Have to find the most-significant bit.
  size_t msb = 0;
  std::cout << "for " << integer;
  // We unroll the loop to reduce complexity.
  if (integer > 0xFFFFFFFF) {
    integer >>= 32;
    msb += 32;
  }
  if (integer > 0xFFFF) {
    integer >>= 16;
    msb += 16;
  }
  if (integer > 0xFF) {
    integer >>= 8;
    msb += 8;
  }
  if (integer > 0xF) {
    integer >>= 4;
    msb += 4;
  }
  if (integer > 0b11) {
    integer >>= 2;
    msb += 2;
  }
  if (integer > 0b1) {
    integer >>= 1;
    msb += 1;
  }
  if (integer > 0) {
    msb += 1;
  }

  // msb uses 1-based counting.
  std::cout  << " msb is " << msb << std::endl;

  for (int i = msb - 1; i >= 0; --i)
    bits->push_back(integer & (1 << i));
}

}  // namespace

std::string Graph::AsGraph6() const {
  // Based on the description of Graph6 here:
  // http://users.cecs.anu.edu.au/~bdm/data/formats.txt
  // Find N(n).
  uint64_t n = vertices_.size();
  std::vector<uint8_t> repr;
  std::vector<bool> N_bits;
  if (n <= 62) {
    repr.push_back(n + 63);
  } else if (n <= 258047) {
    IntegerToBits(n, &N_bits);
    repr.push_back(126);
    Graph6Rx(N_bits, &repr);
  } else if (n <= 68719476735) {
    IntegerToBits(n, &N_bits);
    repr.push_back(126);
    repr.push_back(126);
    Graph6Rx(N_bits, &repr);
  } else {
    std::cerr << "Error: graph too big to represent as Graph6. # nodes = "
              << std::to_string(n) << std::endl;
    return std::string();
  }

  // Assemble the upper triangle of the adjacency matrix.
  // TODO(aryap): Is it better to pre-compute this?
  // Is it better to find every adjacency as needed per node?
  std::vector<std::set<size_t>> adjacent_vertices(vertices_.size());
  for (int i = 0; i < vertices_.size(); ++i) {
    const Vertex *vertex = vertices_.at(i);
    for (const Edge *in_edge : vertex->in()) {
      for (const Vertex *in_vertex : in_edge->in) {
        adjacent_vertices[i].insert(in_vertex->index());
      }
      // TODO(aryap): Are other vertices driven by the same input edge as this
      // one adjacent to this one?
    }
    if (vertex->out() == nullptr) continue;
    for (const Vertex *out_vertex : vertex->out()->out) {
      adjacent_vertices[i].insert(out_vertex->index());
    }
  }

  std::vector<bool> bits;
  for (int i = 0; i < vertices_.size(); ++i) {
    // Skips (0, 0).
    for (int j = 0; j < i; ++j) {
      bits.push_back(adjacent_vertices[j].find(i) != adjacent_vertices[j].end());
    }
  }

  Graph6Rx(bits, &repr);

  return std::string(repr.begin(), repr.end());
}

std::string Graph::AsEdgeListWithWeights() const {
  std::string repr;
  for (const Edge *edge : edges_)
    for (const Vertex *in : edge->in)
      for (const Vertex *out : edge->out)
        repr += std::to_string(in->index()) + " "
                + std::to_string(out->index()) + " "
                + std::to_string(edge->weight) + "\n";
  return repr;
}

} // namespace jsondotrulo
