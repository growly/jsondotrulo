#include "graph.h"

#include <iostream>

#include "edge.h"
#include "vertex.h"

namespace jsondotrulo {

Graph::~Graph() {
  for (const Vertex *vertex : vertices_)
    delete vertex;
  for (const auto &edge : edges_)
    delete edge.second;
}

Edge *Graph::FindOrCreateEdge(const std::string &name) {
  Edge *net;
  auto edge_it = edges_.find(name);
  if (edge_it != edges_.end()) {
    net = edge_it->second;
  } else {
    net = new Edge({name, {}, {}});
    edges_.insert({name, net});
  }
  return net;
}

void Graph::AddVertex(
      const std::vector<std::string> &input_edge_names,
      const std::string &output_edge_name) {
  size_t vertex_index = vertices_.size();
  Vertex *vertex = new Vertex(
      {std::to_string(vertex_index), vertex_index, {}, nullptr, 0.0});
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
        {"IN_" + std::to_string(index), index, {}, nullptr, 0.0});
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
        {"OUT_" + std::to_string(index), index, {}, nullptr, 0.0});
    vertices_.push_back(out_vertex);
    net->out.insert(out_vertex);
    outputs_.push_back(net);
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
}

std::string Graph::AsDOT() const {
  std::string out_str = "digraph " + name_ + "\n{\n";
  for (const auto &iter : edges_) {
    Edge *edge = iter.second;
    for (const Vertex *in : edge->in) {
      for (const Vertex *out : edge->out) {
        out_str += "  " + in->name + " -> " + out->name + " [label=\""
                   + edge->name + "\"]\n";
      }
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

  for (const auto &iter : edges_) {
    Edge *edge = iter.second;
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

} // namespace jsondotrulo
