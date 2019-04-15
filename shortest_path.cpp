/**
 * Visualize shortest paths on a graph using heatmaps.
 */

#include <vector>
#include <fstream>
#include <queue>

#include "utils/sfml_viewer.hpp"
#include "utils/util.hpp"
#include "utils/color.hpp"

#include "graph.hpp"

// define aliases
using GraphType = Graph<int, int>;
using NodeType  = typename GraphType::Node;
using EdgeType  = typename GraphType::Edge;
using NodeID    = typename GraphType::NodeID;
using NodeIter  = typename GraphType::NodeIterator;

// functor comparing node positions relative to reference point (used in nearest_node)
struct CompareNodePositions {
  Point point_;

  bool operator()(const NodeType &node_a, const NodeType &node_b) {
    double dist_a = norm_2(point_ - node_a.position());
    double dist_b = norm_2(point_ - node_b.position());
    return dist_a < dist_b;
  }

  CompareNodePositions(const Point &point) : point_(point) {
  }
};

// find nearest node to a specified point
NodeIter nearest_node(const GraphType& g, const Point& point) {
  return std::min_element(g.node_begin(), g.node_end(), CompareNodePositions(point));
}

 // update graph with shortest path lengths from a root node (Breadth-First Search)
int shortest_path_lengths(GraphType& g, NodeType& root)
{
  // initialize all nodes to -1
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    (*it).value() = -1;
  }
  root.value() = 0; // set root node to 0

  int max_dist = 0; // keeping track of max dist from the root
  std::queue<NodeType> q;
  q.push(root);

  while (!q.empty()) {
    NodeType node = q.front();
    q.pop(); // pop off top of queue

    const int &node_value = node.value();

    // iterate over neighbors
    for (auto ei = node.edge_begin(); ei != node.edge_end(); ++ei) {
      NodeType neighbor = (*ei).node2(); // neighbor node

      int &neighbor_value = neighbor.value();

      // skip neighbor if visited already
      if (neighbor_value > 0) {
        continue;
      }

      // update neighbor value
      neighbor_value = node_value + 1; // constant cost of traversing edge
      q.push(neighbor);

      // keep track of farthest distance seen
      if (neighbor_value > max_dist) {
        max_dist = neighbor_value;
      }
    }
  }

  return max_dist;
}

// functor for generating heat map on graph
struct ColorFunctor {
  int max_value {};

  GraphUtil::Color operator()(const NodeType &node) {
    int node_value = node.value();
    return GraphUtil::Color::make_heat(node_value/(float)max_value);
  }

  ColorFunctor(const int max_value_) : max_value(max_value_) {
  }
};

int main(int argc, char** argv)
{
  // check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // construct graph
  GraphType graph;
  std::vector<GraphType::node_type> nodes;

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
  std::cout << "# Nodes: " << graph.num_nodes() << std::endl;
  std::cout << "# Edges: " << graph.num_edges() << std::endl;

  // launch viewer
  GraphUtil::SFML_Viewer viewer;

  // setup coloring of paths by length
  NodeIter node_iter = nearest_node(graph, Point(0.0, 0.0, 0.0));
  NodeType root = *node_iter;
  int max_length = shortest_path_lengths(graph, root);

  // load graph into viewer
  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(graph.node_begin(), graph.node_end(), ColorFunctor(max_length), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  // center view, enter interactive event loop
  viewer.center_view();
  viewer.event_loop();

  return 0;
}
