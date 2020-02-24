#include "graph.h"

#include <iostream>

#include "edge.h"
#include "vertex.h"

namespace jsondotrulo {

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
  Vertex *vertex =
      new Vertex({std::to_string(vertices_.size()), {}, nullptr, 0.0});
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
    net->in.insert(new Vertex(
          {"EXTERNAL_IN_" + std::to_string(inputs_.size()), {}, nullptr, 0.0}));
    inputs_.push_back(net);
  }
}

void Graph::AddOutputEdges(const std::vector<std::string> &output_edge_names) {
  for (const std::string &output : output_edge_names) {
    Edge *net = FindOrCreateEdge(output);
    net->out.insert(new Vertex(
          {"EXTERNAL_OUT_" + std::to_string(outputs_.size()), {}, nullptr, 0.0}));
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

void Graph::WriteDOT() const {
  std::cout << "digraph " << name_ << std::endl << "{" << std::endl;
  for (const auto &iter : edges_) {
    Edge *edge = iter.second;
    //std::cout << "// " << edge->name
    //          << " in: " << edge->in.size()
    //          << " out: " << edge->out.size() << std::endl;
    for (const Vertex *in : edge->in) {
      for (const Vertex *out : edge->out) {
        std::cout << "  " << in->name << " -> " << out->name
                  << " [label=\"" << edge->name << "\"]" << std::endl;
      }
    }
  }
  std::cout << "}" << std::endl;
}

} // namespace jsondotrulo
