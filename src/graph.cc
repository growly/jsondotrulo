#include <assert.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <set>
#include <sstream>
#include <unordered_map>
#include <utility>

#include <gflags/gflags.h>

#include "edge.h"
#include "graph.h"
#include "path.h"
#include "vertex.h"

namespace jsondotrulo {

DEFINE_bool(show_edge_longest_paths, false,
            "Print the longest path traversing each edge when weighting them.");
DEFINE_string(trace_paths_for, "",
              "Name of vertex whose paths to print.");

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

void Graph::AddVertex(
      const VertexType &type,
      const std::unordered_map<std::string, std::vector<std::string>>
          &input_edges,
      const std::unordered_map<std::string, std::vector<std::string>>
          &output_edges,
      const std::string &instance_of) {
  size_t vertex_index = vertices_.size();
  Vertex *vertex = new Vertex(vertex_index, type);
  vertex->set_instance_of(instance_of);
  vertex->GenerateName();
  for (const auto &input_it : input_edges) {
    const std::string &port = input_it.first;
    for (const auto &input : input_it.second) {
      Edge *in_edge = FindOrCreateEdge(input);
      in_edge->out.insert(vertex);
      vertex->in().push_back(in_edge);
      vertex->in_ports().push_back(port);
    }
  }
  for (const auto &output_it : output_edges) {
    const std::string &port = output_it.first;
    for (const auto &output : output_it.second) {
      Edge *out_edge = FindOrCreateEdge(output);
      out_edge->in.insert(vertex);
      vertex->out().push_back(out_edge);
      vertex->out_ports().push_back(port);
    }
  }
  if (type == VertexType::MODULE) {
    instances_.push_back(vertex);
  }
  vertices_.push_back(vertex);
}

void Graph::ExpandInstances(
    const std::unordered_map<std::string, Graph*> &modules_by_name) {
  while (!instances_.empty()) {
    // Remove expanded instances from vertices_, and then free memory.
    std::set<Edge*> edges_to_delete;

    Vertex *instance = instances_.back();
    instances_.pop_back();
    std::cout << "Expanding instances (" << instances_.size()
              << " remain); current: vertex " << instance->name() << " of type "
              << instance->instance_of() << std::endl;

    // Try to find Graph for instance type in the given map.
    const std::string &instance_of = instance->instance_of();
    auto modules_it = modules_by_name.find(instance_of);
    if (modules_it == modules_by_name.end()) {
      std::cerr << "Warning: Could not find Graph for instance of '"
                << instance_of << "'" << std::endl;
      continue;
    }
    Graph *master = modules_it->second;

    // Maps vertices & nodes in the master Graph to the new copies in the
    // instantiated copy.
    std::unordered_map<Vertex*, Vertex*> vertex_map;
    std::unordered_map<Edge*, Edge*> edge_map;

    // Copy everything.
    for (Edge *old_edge: master->edges_) {
      // TODO(aryap): Generating a 'unique' name for this edge like this
      // doesn't scale. Maybe it's better to just keep number edges from the
      // last edge in the top module. So this name should be unique. Otherwise
      // we'd just have to reproduce some of the FindOrCreateEdge functionality
      // and not rely on unique names. But we'd still need a unique name. So
      // whatever.
      std::string new_edge_name = instance->name() + "." + old_edge->name;
      Edge *new_edge = FindOrCreateEdge(new_edge_name);
      edge_map[old_edge] = new_edge;
      //std::cout << new_edge_name << " created " << new_edge << std::endl;
    }

    for (Vertex *old_vertex : master->vertices_) {
      // Skip copying input/output pins.
      if (old_vertex->type() == VertexType::IN_PIN ||
          old_vertex->type() == VertexType::OUT_PIN)
        continue;
      // Get new index for copied instance.
      size_t new_index = vertices_.size();
      Vertex *new_vertex = new Vertex(new_index, old_vertex->type());
      new_vertex->set_instance_of(old_vertex->instance_of());
      new_vertex->GenerateName();
      if (new_vertex->type() == VertexType::MODULE)
        instances_.push_back(new_vertex);
      new_vertex->out_ports().insert(new_vertex->out_ports().begin(),
                                     old_vertex->out_ports().begin(),
                                     old_vertex->out_ports().end());
      for (Edge *old_out_edge : old_vertex->out()) {
        auto edge_map_it = edge_map.find(old_out_edge);
        // There should be a 1:1 correspondence of old and new edges!
        assert(edge_map_it != edge_map.end());
        Edge *new_edge = edge_map_it->second;
        new_edge->in.insert(new_vertex);
        new_vertex->out().push_back(new_edge);
      }
      new_vertex->in_ports().insert(new_vertex->in_ports().begin(),
                                    old_vertex->in_ports().begin(),
                                    old_vertex->in_ports().end());
      for (Edge *old_in_edge : old_vertex->in()) {
        auto edge_map_it = edge_map.find(old_in_edge);
        // There should be a 1:1 correspondence of old and new edges!
        assert(edge_map_it != edge_map.end());
        Edge *new_edge = edge_map_it->second;
        //std::cout << "v " << new_vertex << " (" << new_vertex->name() << ") added to " << new_edge << " " << new_edge->name << std::endl;
        new_edge->out.insert(new_vertex);
        new_vertex->in().push_back(new_edge);
      }
      vertices_.push_back(new_vertex);
    }

    // To avoid creating and maintaining this structure for every vertex, we
    // create it here.
    std::unordered_map<std::string, std::vector<Edge*>> instance_input_ports;
    for (int i = 0; i < instance->in().size(); ++i) {
      const std::string &port_name = instance->in_ports().at(i);
      Edge *edge = instance->in().at(i);
      instance_input_ports[port_name].push_back(edge);
    }

    // Now map inputs to the old vertex to those in the new vertex.
    for (auto &pair : instance_input_ports) {
      const std::string &port_name = pair.first;
      std::vector<Edge*> &instance_inputs = pair.second;

      // Find corresponding edges in master Graph for this instance.
      auto master_inputs_it = master->inputs_.find(port_name);
      if (master_inputs_it == master->inputs_.end()) {
        std::cerr << "Error! Port " << port_name << " of instance vertex "
                  << instance->name() << " was not found in master Graph for "
                  << instance->instance_of() << std::endl;
        continue;
      }
      std::vector<Edge*> &master_inputs = master_inputs_it->second;

      for (int i = 0; i < instance_inputs.size(); ++i) {
        Edge *external = instance_inputs[i];
        Edge *master = master_inputs[i];
        auto internal_edge_it = edge_map.find(master);
        assert(internal_edge_it != edge_map.end());
        Edge *internal = internal_edge_it->second;

        // Remove the instance vertex (the one of type MODULE that is to
        // expanded into the full module); replace it with the internal vertex
        // within the graph to which that edge would connect.
        //std::cout << "internal e " << internal << std::endl;
        //std::cout << "[inputs] vertex" <<  instance << " (" << instance->name() << ") being removed from edge "
        //          << external << " (" << external->name << ") " << std::endl;

        size_t erased = external->out.erase(instance);
        // Since external->out is a set but we track multiple connections to
        // the same net through a vector in the input list, this may be
        // redundant on successive iterations.
   
        // Merge the input/output sets between the internal/external edge.
        external->in.insert(internal->in.begin(), internal->in.end());
        external->out.insert(internal->out.begin(), internal->out.end());
        
        // Replace every occurrence of the temporary internal edge - at the
        // vertices which are both inputs and outputs to the edge - with the
        // external edge subsuming it.
        for (Vertex *vertex : internal->in) {
          // There might be multiple connections to this edge!
          //std::cout << "[inputs] replacing e " << internal << "  with e " << external
          //          << " on v " << vertex << std::endl;
          std::replace(vertex->in().begin(), vertex->in().end(), internal, external);
        }
        for (Vertex *vertex : internal->out) {
          // There might be multiple connections to this edge!
          //std::cout << "[inputs] replacing e " << internal << "  with e " << external
          //          << " on v " << vertex << std::endl;
          std::replace(vertex->in().begin(), vertex->in().end(), internal, external);
        }
        internal->in.clear();
        internal->out.clear();
        edges_to_delete.insert(internal);
      }
    }

    // TODO(aryap): Refactor this and the above stanza into a separate
    // function.
    //
    // And then map outputs from the old vertex to those in the elaborated
    // module.
    std::unordered_map<std::string, std::vector<Edge*>> instance_output_ports;
    for (int i = 0; i < instance->out().size(); ++i) {
      const std::string &port_name = instance->out_ports().at(i);
      Edge *edge = instance->out().at(i);
      instance_output_ports[port_name].push_back(edge);
    }

    // Now map outputs from the old vertex to those in the new vertex.
    for (auto &pair : instance_output_ports) {
      const std::string &port_name = pair.first;
      std::vector<Edge*> &instance_outputs = pair.second;

      // Find corresponding edges in master Graph for this instance.
      auto master_outputs_it = master->outputs_.find(port_name);
      if (master_outputs_it == master->outputs_.end()) {
        std::cerr << "Error! Port " << port_name << " of instance vertex "
                  << instance->name() << " was not found in master Graph for "
                  << instance->instance_of() << std::endl;
        continue;
      }
      std::vector<Edge*> &master_outputs = master_outputs_it->second;

      for (int i = 0; i < instance_outputs.size(); ++i) {
        Edge *external = instance_outputs[i];
        Edge *master = master_outputs[i];
        auto internal_edge_it = edge_map.find(master);
        assert(internal_edge_it != edge_map.end());
        Edge *internal = internal_edge_it->second;

        // Remove the instance vertex (the one of type MODULE that is to
        // expanded into the full module); replace it with the internal vertex
        // within the graph to which that edge would connect.
        //std::cout << "internal e " << internal << std::endl;
        //std::cout << "[outputs] vertex" <<  instance << " (" << instance->name() << ") being removed from edge "
        //          << external << " (" << external->name << ") " << std::endl;
        size_t erased = external->in.erase(instance);

        // Merge the input/output sets between the internal/external edge.
        external->in.insert(internal->in.begin(), internal->in.end());
        external->out.insert(internal->out.begin(), internal->out.end());
        
        
        for (Vertex *vertex : internal->in) {
        //  std::cout << "[outputs] replacing e " << internal << "  with e " << external
        //            << " on v " << vertex << std::endl;
          std::replace(vertex->out().begin(), vertex->out().end(), internal, external);
        }
        for (Vertex *vertex : internal->out) {
        //  std::cout << "[outputs] replacing e " << internal << "  with e " << external
        //            << " on v " << vertex << std::endl;
          std::replace(vertex->in().begin(), vertex->in().end(), internal, external);
        }
        internal->in.clear();
        internal->out.clear();
        edges_to_delete.insert(internal);
      }
    }

    // Remove the instance vertex and all the internal input/output edges we
    // copied.
    auto vertex_iter = std::find(
        vertices_.begin(), vertices_.end(), instance);
    assert(vertex_iter != vertices_.end());
    vertices_.erase(vertex_iter);

    assert(std::find(
        instances_.begin(), instances_.end(), instance) == instances_.end());
    
    // Have to reassign indices.
    for (size_t i = 0; i < vertices_.size(); ++i) {
      vertices_[i]->set_index(i);
      // TODO(aryap): Should Vertices be renamed?
      //vertices_[i]->GenerateName();
    }

    for (Edge *edge : edges_to_delete) {
      //std::cout << "deleting e " << edge << " " << edge->name << std::endl;
      assert(edge->in.empty());
      assert(edge->out.empty());

      auto edges_by_name_it = edges_by_name_.find(edge->name);
      assert(edges_by_name_it != edges_by_name_.end());
      edges_by_name_.erase(edges_by_name_it);

      auto edge_it = std::find(edges_.begin(), edges_.end(), edge);
      assert(edge_it != edges_.end());
      edges_.erase(edge_it);

      // This is optimistic.
      auto edges_by_partition_it = edges_by_partition_.find(edge->partition);
      assert(edges_by_partition_it != edges_by_partition_.end());

      std::set<Edge*> &edge_set = edges_by_partition_it->second;
      auto edge_set_it = edge_set.find(edge);
      assert(edge_set_it != edge_set.end());
      edge_set.erase(edge_set_it);

      delete edge;
    }

    assert(edges_by_name_.size() == edges_.size());

    // Reassign Edge indices.
    for (size_t i = 0; i < edges_.size(); ++i) {
      edges_[i]->index = i;
      // Edge names are NOT modified.
    }
  }
}

void Graph::UpdateEdgeWeightsForPath(
    Path *path,
    std::set<Path*> *critical_paths,
    std::unordered_map<Edge*, Path*> *critical_path_by_edge,
    bool *true_on_update) {
  double path_cost = path->Cost();
  for (auto &pair : *path) {
    Edge *edge = pair.second;
    if (edge == nullptr) continue;
    if (path_cost >= edge->weight) {
      edge->weight = path_cost;
      (*critical_path_by_edge)[edge] = path;
      critical_paths->insert(path);
      if (true_on_update != nullptr) *true_on_update = true;
    }
  }
}

// Traverse the tree. At every node, record the longest path from that point.
// This avoids repeating searches.
void Graph::WeightCombinatorialPaths2() {
  std::set<Path*> paths;
  std::unordered_map<Vertex*, Path*> descendant_by_vertex;
  std::unordered_map<Edge*, Path*> descendant_by_edge;

  std::set<Path*> critical_paths;
  std::unordered_map<Edge*, Path*> critical_path_by_edge;

  std::vector<std::pair<Edge*, Path*>> to_visit;
  
  size_t num_found = 0;

  // Seed start paths.
  for (Vertex *start_vertex : vertices_) {
    if (!start_vertex->IsSynchronous()) continue;

    for (Edge *start_edge : start_vertex->out()) {
      Path *start_path = new Path{{start_vertex, start_edge}};
      paths.insert(start_path);

      to_visit.push_back(std::make_pair(start_edge, start_path));
    }
  }

  while (!to_visit.empty()) {
    Edge *current = to_visit.back().first;
    Path *path = to_visit.back().second;
    to_visit.pop_back();

    if (path->front().first->name() == FLAGS_trace_paths_for) {
      path->Print();
    }

    auto descendant_it = descendant_by_edge.find(current);
    if (descendant_it != descendant_by_edge.end()) {
      // We already have the longest descendant for this edge, so we can
      // assemble the longest path through it.
      Path *final_path = new Path(*path);
      final_path->Append(*descendant_it->second);

      // Update edge weights.
      bool defer_delete = false;
      UpdateEdgeWeightsForPath(final_path,
                               &critical_paths,
                               &critical_path_by_edge,
                               &defer_delete);
      if (!defer_delete) delete final_path;

      // We're done with this edge now; move on.
      std::cout << "continuing: edge " << current->name << " has known descendants w: "
                << descendant_it->second->Cost() << "; "
                << descendant_it->second->AsString() << std::endl;
      continue;
    }

    bool current_edge_is_terminal = true;

    double vertex_max_cost = -1.0;
    Vertex *max_cost_vertex = nullptr;
    // If this is non-null, it indicates that the current edge is terminal
    // because the max_cost_vertex has known descendants. Then it points to the
    // longest path of which max_cost_vertex is the root.
    Path *max_cost_vertex_path = nullptr;

    for (Vertex *next_vertex : current->out) {
      // Ignore loops.
      if (path->ContainsVertex(next_vertex))
        continue;

      // If the cost of the subgraph below this vertex is already known, use
      // it.
      auto descendant_by_vertex_it = descendant_by_vertex.find(next_vertex);
      if (descendant_by_vertex_it != descendant_by_vertex.end()) {
        Path *terminal_path = descendant_by_vertex_it->second;

        if (terminal_path->Cost() > vertex_max_cost) {
          max_cost_vertex_path = terminal_path;
          max_cost_vertex = next_vertex;
          vertex_max_cost = terminal_path->Cost();
        }

        std::cout << "vertex "  << next_vertex->name()
                  << " has known descendant path w: "
                  << terminal_path->Cost()
                  << " " << terminal_path->AsString() << std::endl;

        Path *final_path = new Path(*path);
        final_path->Append(*terminal_path);

        bool defer_delete = false;
        UpdateEdgeWeightsForPath(final_path,
                                 &critical_paths,
                                 &critical_path_by_edge,
                                 &defer_delete);
        if (!defer_delete) delete final_path;
        continue;
      }

      if (next_vertex->Cost() > vertex_max_cost) {
        max_cost_vertex_path = nullptr;
        max_cost_vertex = next_vertex;
        vertex_max_cost = next_vertex->Cost();
      }

      // Paths end at flip-flops and latches and such.
      if (next_vertex->IsSynchronous()) {
        bool defer_delete = false;
        Path *final_path = new Path(*path);
        assert(!final_path->ContainsVertex(next_vertex));
        final_path->Append(
            std::make_pair(next_vertex, nullptr));

        // It takes too much memory to save these, so we process them now and
        // move on. Unless we have to.
        ++num_found;

        UpdateEdgeWeightsForPath(final_path,
                                 &critical_paths,
                                 &critical_path_by_edge,
                                 &defer_delete);

        // Unless we're doing some intense admin work, we delete the path
        // now that we're done with it.
        if (!FLAGS_show_edge_longest_paths || !defer_delete)
          delete final_path;
        continue;
      }

      current_edge_is_terminal = false;

      // If the path-so-far isn't ending, extend a copy of the path with
      // the next edge to follow.
      //
      // Also check if this vertex is terminal, which happens when all of its
      // child edges are terminal.
      bool next_vertex_is_terminal = true;
      double edge_max_cost = -1.0;
      Path *max_cost_edge_path = nullptr;
      Edge *max_cost_edge = nullptr;
      for (Edge *next_edge : next_vertex->out()) {
        descendant_it = descendant_by_edge.find(next_edge);
        if (descendant_it == descendant_by_edge.end()) {
          next_vertex_is_terminal = false;
        } else {
          Path *descendant_path = descendant_it->second;
          if (descendant_path->Cost() > edge_max_cost) {
            max_cost_edge_path = descendant_path;
            max_cost_edge = next_edge;
            edge_max_cost = descendant_path->Cost();
          }
        }

        // Don't add a duplicate edge to the path.
        if (path->ContainsEdge(next_edge))
          continue;
        // Create a new path for the each next-hop.
        Path *new_path = new Path(*path);

        new_path->Append(std::make_pair(next_vertex, next_edge));
        paths.insert(new_path);
        to_visit.push_back(std::make_pair(next_edge, new_path));

        std::cout << "expanding search to: ";
        new_path->Print();
      }

      if (next_vertex_is_terminal && max_cost_edge_path != nullptr) {
        std::cout << "vertex " << next_vertex->name()
            << " has all-known descendants with longst path ";
        Path *terminal_path = new Path({next_vertex, max_cost_edge});
        terminal_path->Append(*max_cost_edge_path);
        terminal_path->Print();
        descendant_by_vertex[next_vertex] = terminal_path;
      }
    }

    if (current->out.empty()) {
      std::cout << "edge " << current->name << " has no out vertices" << std::endl;
    }

    if (current_edge_is_terminal && !current->out.empty()) {
      // All of the vertices out of this edge are terminal, so record the
      // longest descendant. At this point we need to find the highest-cost
      // final node and use that as the terminal node.
      Path *terminal_path = max_cost_vertex_path == nullptr ?
          new Path(std::pair<Vertex*, Edge*>{max_cost_vertex, nullptr}) :
          new Path(*max_cost_vertex_path);
      descendant_by_edge[current] = terminal_path;
      std::cout << "edge " << current->name << " is terminal: ";
      terminal_path->Print();
    }

    delete path;
  }

  if (FLAGS_show_edge_longest_paths) {
    for (Edge *edge : edges_) {
      std::cout << "(" << edge->name << " wt: " << edge->weight << ") [";
      auto it = critical_path_by_edge.find(edge);
      if (it == critical_path_by_edge.end()) {
        std::cout << "none]" << std::endl;
        continue;
      }
      Path *path = it->second;
      path->Print();
    }

    for (Path *path : critical_paths)
      delete path;
  }

  // TODO(aryap): Delete descendant paths.
  for (const auto &pair : descendant_by_edge) {
    std::cout << "longest descendant of edge " << pair.first->name
              << " w: " << pair.second->Cost() << " "
              << pair.second->AsString() << std::endl;
  }

  for (const auto &pair : descendant_by_vertex) {
    std::cout << "longest descendant of vertex " << pair.first->name()
              << " w: " << pair.second->Cost() << " "
              << pair.second->AsString() << std::endl;
  }

}

// TODO(aryap): I'm an idiot.  This is stupid. We don't need to enumerate every
// possible path between every two endpoints. We need to find the the longest
// path between any two points once, which is why we need to use Floyd-Warshall
// (again). Dijkstra/Floyd-Warshall would stop us repeating searching on
// subgraphs that we'd already traversed, that's why this is stupid. I was
// doing that before without considering longer alternatives.
//
//
void Graph::WeightCombinatorialPaths() {
  // The goal is to find every path from one flip flop to another flip flop,
  // through whatever other vertices.

  std::cout << "Finding paths between synchronous elements" << std::endl;

  // Paths must be deleted at the end of this scope!
  std::set<Path*> paths;
  size_t num_found = 0;
  size_t i = 0;
  size_t update_interval =
      std::max(vertices_.size() / 100, static_cast<size_t>(1));

  std::set<Path*> critical_paths;
  std::unordered_map<Edge*, Path*> critical_path_by_edge;

  for (Vertex *start_vertex : vertices_) {
    ++i;
    if (i % update_interval == 0) {
      std::cout << "Finding paths: " << i << "/" << vertices_.size()
                << " start points checked; "
                << num_found << " paths found" << std::endl;
    }

    if (!start_vertex->IsSynchronous()) continue;

    for (Edge *start_edge : start_vertex->out()) {
      Path *start_path = new Path{{start_vertex, start_edge}};
      paths.insert(start_path);

      std::vector<std::pair<Edge*, Path*>> to_visit;
      to_visit.push_back(std::make_pair(start_edge, start_path));

      // We do not care for the other drivers of the same net. Presumably they
      // are not on the path.
      while (!to_visit.empty()) {
        Edge *current = to_visit.back().first;
        Path *path = to_visit.back().second;
        to_visit.pop_back();

        bool defer_delete = false;

        for (Vertex *next_vertex : current->out) {
          if (path->ContainsVertex(next_vertex))
            continue;

          // Paths end at flip-flops.
          if (next_vertex->IsSynchronous()) {
            Path *final_path = new Path(*path);
            assert(!final_path->ContainsVertex(next_vertex));
            final_path->Append(
                std::make_pair(next_vertex, nullptr));

            // It takes too much memory to save these, so we process them now
            // and move on. Unless we have to.
            ++num_found;
            double path_cost = final_path->Cost();

            for (auto &pair : *final_path) {
              Edge *edge = pair.second;
              if (edge == nullptr) continue;

              if (path_cost >= edge->weight) {
                edge->weight = path_cost;
                critical_path_by_edge[edge] = final_path;
                critical_paths.insert(final_path);
                defer_delete = true;
              }
            }

            // Unless we're doing some intense admin work, we delete the path
            // now that we're done with it.
            if (!FLAGS_show_edge_longest_paths || !defer_delete)
              delete final_path;
            continue;
          }

          // If the path-so-far isn't ending, extend a copy of the path with
          // the next edge to follow.
          for (Edge *next_edge : next_vertex->out()) {
            // Don't add a duplicate edge to the path.
            if (path->ContainsEdge(next_edge))
              continue;
            // Create a new path for the each next-hop.
            Path *new_path = new Path(*path);

            new_path->Append(std::make_pair(next_vertex, next_edge));
            paths.insert(new_path);
            to_visit.push_back(std::make_pair(next_edge, new_path));

            if (start_vertex->name() == FLAGS_trace_paths_for)
              new_path->Print();
          }
        }

        delete path;
      }
    }
  }

  std::cout << "Finding paths: done" << std::endl;

  if (FLAGS_show_edge_longest_paths) {
    for (Edge *edge : edges_) {
      std::cout << "(" << edge->name << " wt: " << edge->weight << ") [";
      auto it = critical_path_by_edge.find(edge);
      if (it == critical_path_by_edge.end()) {
        std::cout << "none]" << std::endl;
        continue;
      }
      Path *path = it->second;
      path->Print();
    }

    for (Path *path : critical_paths)
      delete path;
  }
}

void Graph::AddInputEdges(
    const std::unordered_map<std::string, std::vector<std::string>>
        &input_edges) {
  for (const auto &pair : input_edges) {
    const std::string &port_name = pair.first;
    for (const std::string &input : pair.second) {
      Edge *net = FindOrCreateEdge(input);
      size_t index = vertices_.size();
      Vertex *in_vertex = new Vertex(index, VertexType::IN_PIN);
      in_vertex->GenerateName();
      vertices_.push_back(in_vertex);
      net->in.insert(in_vertex);
      inputs_[port_name].push_back(net);
    }
  }
}

void Graph::AddOutputEdges(
    const std::unordered_map<std::string, std::vector<std::string>>
        &output_edges) {
  for (const auto &pair : output_edges) {
    const std::string &port_name = pair.first;
    for (const std::string &output : pair.second) {
      Edge *net = FindOrCreateEdge(output);
      size_t index = vertices_.size();
      Vertex *out_vertex = new Vertex(index, VertexType::OUT_PIN);
      out_vertex->GenerateName();
      vertices_.push_back(out_vertex);
      net->out.insert(out_vertex);
      outputs_[port_name].push_back(net);
    }
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
  for (const auto &pair : inputs_) {
    std::cout << "\t\t" << pair.first << std::endl;
    for (const Edge *edge : pair.second)
      std::cout << "\t\t\t" << edge->name << std::endl;
  }
  std::cout << "\toutputs: " << outputs_.size() << std::endl;
  for (const auto &pair : outputs_) {
    std::cout << "\t\t" << pair.first << std::endl;
    for (const Edge *edge : pair.second)
      std::cout << "\t\t\t" << edge->name << std::endl;
  }

  std::cout << "\tvertices: " << vertices_.size() << std::endl;
  std::cout << "\tof which instances: " << instances_.size() << std::endl;
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
    for (const Edge *out_edge : vertex->out()) {
      for (const Vertex *out_vertex : out_edge->out) {
        adjacent_vertices[i].insert(out_vertex->index());
      }
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
      for (const Vertex *out : edge->out) {
        std::string line = in->name() + " " + out->name() + " "
                           + std::to_string(edge->weight) + "\n";
        std::cout << "edge: " << edge->name << " w: " << line;
        repr += line;
      }
  return repr;
}

} // namespace jsondotrulo
