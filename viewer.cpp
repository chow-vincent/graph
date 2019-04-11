 /**
  * Visualize graphs, with nodes represented as 3D points in space.
  */

#include <fstream>

#include "utils/sfml_viewer.hpp"
#include "utils/util.hpp"

#include "graph.hpp"

int main(int argc, char** argv)
{
  // check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // define aliases
  using GraphType = Graph<double, double>;
  using NodeType  = typename GraphType::node_type;

  // construct graph
  GraphType graph;
  std::vector<NodeType> nodes;

  // nodes_file from first input arg
  std::ifstream nodes_file(argv[1]);

  // each line of node_file is a 3D point, add to graph
  Point p;
  while (GraphUtil::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // tets_file from second input arg
  std::ifstream tets_file(argv[2]);

  // each line of tets_file as four ints referring to node IDs
  std::array<int,4> t;
  while (GraphUtil::getline_parsed(tets_file, t))
    for (unsigned i = 1; i < t.size(); ++i)
      for (unsigned j = 0; j < i; ++j)
        graph.add_edge(nodes[t[i]], nodes[t[j]]);

  // print num of nodes and edges
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // launch viewer
  GraphUtil::SFML_Viewer viewer;

  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  // center view, enter interactive event loop
  viewer.center_view();
  viewer.event_loop();

  return 0;
}
