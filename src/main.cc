#include <assert.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <unordered_map>

#include "c_make_header.h"

#include "graph.h"
#include "vertex.h"
#include "edge.h"

#include <gflags/gflags.h>
#include <nlohmann/json.hpp>

DEFINE_string(dot_file, "",
              "Path prefix to DOT-format output file. If empty, no DOT output "
              "will be provided.");
DEFINE_string(m_file, "",
              "Path prefix to MATLAB m-file output. If empty, no MATLAB output "
              "will be provided.");
DEFINE_string(hmetis_file, "",
              "Path prefix to hMETIS input file to write. If empty, no MATLAB"
              " output will be provided.");
DEFINE_string(graph6_file, "",
              "Path to prefix to Graph6 file to write. If empty, no Graph6 "
              " output will be provided.");
DEFINE_string(edge_list, "",
              "Path to prefix for edge-list file to write. If empty, no"
              " edge-list output will be provided.");
DEFINE_bool(print, false, "Print details about the module graph.");
DEFINE_bool(expand_instances, true, "Expand instance definitions.");

// TODO(aryap): The ordering has to be deterministic between runs if we're just
// going to ingest this metadata like this:
DEFINE_string(hmetis_partition_file, "",
              "Path to hMETIS-output partition file, containing a partition "
              "number per line, indicating the allocating partition number "
              "for the hyperedge of the same number");

DEFINE_string(top, "", "Name of the top verilog module.");

void WriteFile(const std::string &file_name, const std::string &content) {
  std::ofstream out_file;
  out_file.open(file_name, std::ios::out | std::ios::trunc);
  out_file << content;
  out_file.close();
}

std::string ReadFile(const std::string &file_name) {
  std::ifstream in_file;
  in_file.open(file_name, std::ios::in);
  // This double-parantheses crap circumvents the "most vexing parse" problem.
  std::string content((std::istreambuf_iterator<char>(in_file)),
                       std::istreambuf_iterator<char>());
  in_file.close();
  return content;
}

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cout << argv[0] << " Version " <<  jsondotrulo_VERSION_MAJOR << "."
              << jsondotrulo_VERSION_MINOR << std::endl;
  }
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  // Read json file.
  std::ifstream i(argv[1]);
  nlohmann::json j;
  i >> j;

  std::unordered_map<std::string, jsondotrulo::Graph*> modules_by_name;


  // Construct Graphs from ingested JSON.
  if (j.find("modules") != j.end()) {
    for (const auto &module_json : j["modules"].items()) {
      jsondotrulo::Graph *g = new jsondotrulo::Graph(module_json.key());
      // Cells become vertices. Each has a list of ports. Each port has an
      // entry in a separate object that determines whether it is input or
      // output, which we have to consult. Ports are connected to multiple
      // nets, denoted by integers. We also need to check the "type:" field, to
      // see if this is a LUT or a flip-flop.
      if (module_json.value().find("cells") != module_json.value().end()) {
        for (const auto &cells_json_it : module_json.value()["cells"].items()) {
          std::unordered_map<std::string, std::vector<std::string>> in_ports;
          std::unordered_map<std::string, std::vector<std::string>> out_ports;
          const auto &cells_json = cells_json_it.value();
          // Have to find connections with type "input" or "output":
          for (const auto &port_json : cells_json["connections"].items()) {
            const std::string &port_name = port_json.key();
            if (cells_json.find("port_directions") == cells_json.end())
              continue;
            const auto &dir_json = cells_json["port_directions"];
            if (dir_json.find(port_name) == dir_json.end())
              continue;
            if (dir_json[port_name] == "input") {
              for (const auto &net_json : port_json.value()) {
                in_ports[port_name].push_back(net_json.dump());
              }
            } else if (dir_json[port_name] == "output") {
              for (const auto &net_json : port_json.value()) {
                out_ports[port_name].push_back(net_json.dump());
              }
            }
          }
          // assert(out_ports.size() == 1); // TODO(aryap): make sure this still works.
          // Figure out what kind of node we have.
          jsondotrulo::VertexType type = jsondotrulo::VertexType::UNKNOWN;
          const std::string &cell_type = cells_json["type"];
          if (cell_type[0] != '$') {
            type = jsondotrulo::VertexType::MODULE;
          } else if (cell_type.find("$lut") != std::string::npos) {
            type = jsondotrulo::VertexType::LUT;
          } else if (cell_type.find("DFF") != std::string::npos) {
            type = jsondotrulo::VertexType::FLIP_FLOP;
          } else {
            std::cout << "Warning, could not identify cell type: " << cell_type
                      << std::endl;
          }
            
          g->AddVertex(type, in_ports, out_ports, cell_type);
        }
      }
      // Find input and output ports. Each port's name is a key, and the object
      // for that key contains both the direction of the port and the nets to
      // which it maps. So we have to look up whether the port is input or
      // output and then iterate over the named nets.
      if (module_json.value().find("ports") != module_json.value().end()) {
        std::unordered_map<std::string, std::vector<std::string>> in_ports;
        std::unordered_map<std::string, std::vector<std::string>> out_ports;
        for (const auto &ports_json_it : module_json.value()["ports"].items()) {
          const std::string &port_name = ports_json_it.key();
          if (ports_json_it.value().find(
                  "direction") == ports_json_it.value().end()) continue;
          const std::string &direction = ports_json_it.value()["direction"];
          if (ports_json_it.value().find(
                  "bits") == ports_json_it.value().end()) continue;
          if (direction == "input") {
            for (const auto &bits_json : ports_json_it.value()[
                "bits"].items()) {
              in_ports[port_name].push_back(bits_json.value().dump());
            }
          } else if (direction == "output") {
            for (const auto &bits_json : ports_json_it.value()[
                "bits"].items()) {
              // This used to be nice. I used to be able to assume every value
              // was an int and use bits_json.value().get<int>(). But SOME
              // values are anot ints, they are strings. FFS.
              out_ports[port_name].push_back(bits_json.value().dump());
            }
          }
        }
        g->AddInputEdges(in_ports);
        g->AddOutputEdges(out_ports);
      }
      if (!FLAGS_hmetis_partition_file.empty()) {
        std::cout << g->name() << ": adding partition data from hMETIS file "
                  << FLAGS_hmetis_partition_file << std::endl;
        g->ReadHMETISPartitions(ReadFile(FLAGS_hmetis_partition_file));
      }
      modules_by_name.insert({g->name(), g});
      if (FLAGS_print) {
        g->Print();
      }
    }
  }
  if (modules_by_name.empty()) {
    std::cerr << "No modules were ingested!" << std::endl;
    return EXIT_FAILURE;
  } else {
    std::cout << "Ingested " << modules_by_name.size()
              << " module description(s)." << std::endl;
  }

  jsondotrulo::Graph *top;
  // After initial module read, we can instantiate child modules. But first we
  // have to determine our starting point (top).
  if (FLAGS_top.empty()) {
    // TODO(aryap): Or maybe make a guess?
    std::cerr << "No top module specified." << std::endl;
    return EXIT_FAILURE;
  }
  auto module_it = modules_by_name.find(FLAGS_top);
  if (module_it != modules_by_name.end()) {
    top = module_it->second;
  } else {
    std::cerr << "Top module \"" << FLAGS_top << "\" not found." << std::endl;
    return EXIT_FAILURE;
  }
  if (FLAGS_expand_instances) {
    top->ExpandInstances(modules_by_name);
    if (FLAGS_print) {
      top->Print();
    }
  }
  top->WeightCombinatorialPaths();
  // Graph object should now be complete. Write different formats now.
  // Except if we weren't given a --top argument, we probably want to avoid
  // the complications of writing the module name to the filename:
  std::string top_part = FLAGS_top.empty() ? "" : "." + top->name();
  if (!FLAGS_dot_file.empty()) {
    std::string file_name = FLAGS_dot_file + top_part + ".gv";
    std::cout << top->name() << ": wrote " << file_name << std::endl;
    WriteFile(file_name, top->AsDOT());
  }
  if (!FLAGS_m_file.empty()) {
    std::string file_name = FLAGS_m_file + top_part + ".m";
    std::cout << top->name() << ": wrote " << file_name << std::endl;
    WriteFile(file_name, top->AsMFile());
  }
  if (!FLAGS_hmetis_file.empty()) {
    std::string file_name = FLAGS_hmetis_file + top_part + ".hmetis";
    std::cout << top->name() << ": wrote " << file_name << std::endl;
    WriteFile(file_name, top->AsHMETIS());
  }
  if (!FLAGS_graph6_file.empty()) {
    std::string file_name = FLAGS_graph6_file + top_part + ".g6";
    std::cout << top->name() << ": wrote " << file_name << std::endl;
    WriteFile(file_name, top->AsGraph6());
  }
  if (!FLAGS_edge_list.empty()) {
    std::string file_name = FLAGS_edge_list + top_part + ".edges";
    std::cout << top->name() << ": wrote " << file_name << std::endl;
    WriteFile(file_name, top->AsEdgeListWithWeights());
  }

  return EXIT_SUCCESS;
}
